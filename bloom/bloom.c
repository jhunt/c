#include <stdint.h>
#include <string.h>
#include <math.h>
#include "insist.h"

static uint32_t s_djb(const char *str, size_t len)
{
	unsigned int hash = 5381;
	size_t i;

	for (i = 0; i < len; i++)
		hash = ((hash << 5) + hash) + str[i]; /* h * 33 + c */

	return hash;
}

static const uint32_t C1 = 0xcc9e2d51;
static const uint32_t C2 = 0x1b873593;
static const uint32_t R1 = 15;
static const uint32_t R2 = 13;
static const uint32_t M  = 5;
static const uint32_t N  = 0xe6546b64;

/* rotate (x) to the left (y) bits */
#define ROL(x,y) ((x) << (y) | ((x) >> (32 - (y))))

static uint32_t s_murmur3(const char *str, size_t len, uint32_t seed)
{
	uint32_t hash, k;
	unsigned long nhead;
	const uint32_t *head;
	const uint8_t  *tail;
	int i;

	hash = seed;
	nhead = len / 4;
	head = (const uint32_t *)str;

	for (i = 0; i < nhead; i += 4) {
		k = ROL(head[i] * C1, R1) * C2;
		hash = ROL(hash^k, R2) * M + N;
	}

	k = 0;
	tail = (uint8_t *)str + (len / 4);
	switch (len & 3) {
	case 3: k ^= tail[2] << 16;
	case 2: k ^= tail[1] << 8;
	case 1: k ^= tail[0];
	        hash ^= ROL(k * C1, R1) * C2;
	}

	hash ^= (uint32_t)len;
	hash ^= (hash >> 16);
	hash *= 0x85ebca6b;
	hash ^= (hash >> 13);
	hash *= 0xc2b2ae35;
	hash ^= (hash >> 16);

	return hash;
}

#define WORDSIZE sizeof(int) * 8

struct _bloom_t {
	unsigned int m;  /* how wide is the filter, in bits?    */
	unsigned int mn; /* m:n ratio, used for fp calcualtions */
	unsigned int k;  /* how many times do we hash each key? */
	int b[]; /* horribly inefficient */
};
typedef struct _bloom_t *bloom_t;

/* Create a new Bloom Filter, of n bits.  m and k will be chosen
   appropriate to the desired false positive error rate. */
bloom_t bloom_new(unsigned int n, unsigned int mn)
{
	insist(n  > 0, "bloom_new(n, mn) is only defined for n > 0");
	insist(mn > 1, "bloom_new(n, mn) is only defined for mn > 1");
	mn *= 8;

	bloom_t bloom = calloc(1, sizeof(bloom_t) + sizeof(int) * n / WORDSIZE);
	insist(bloom != NULL, "memory allocation failed");

	bloom->m  = n;
	bloom->mn = mn;
	bloom->k  = (unsigned int)(round(log(2*mn) / log(2)));
	return bloom;
}

void free_bloom(bloom_t bloom)
{
	free(bloom);
}

double bloom_fp(bloom_t bloom)
{
	insist(bloom != NULL, "bloom_fp(NULL) is undefined");
	return pow(1 - exp(-1.0 * bloom->k / bloom->mn), bloom->k);
}

void bloom_diag(bloom_t bloom, FILE *io, const char *prefix)
{
	int i;
	if (prefix == NULL) prefix = "";

	fprintf(io, "%s[bloom %p]\n", prefix, bloom);
	fprintf(io, "%s m = %u, k = %u, e = %lf\n", prefix, bloom->m, bloom->k, bloom_fp(bloom));

	if (bloom->m <= (64 * 64)) {
		int wrap = 64;
		for (i = 0; i < bloom->m; i++) {
			if (i % wrap == 0) fprintf(io, "%s%s  [ ", (i ? "]\n" : ""), prefix);
			fprintf(io, "%u ", bloom->b[i]);
		}
		fprintf(io, "]\n");
	}
}

void bloom_set(bloom_t bloom, const char *key, size_t len)
{
	insist(bloom != NULL, "bloom_set() is undefined for NULL bloom filters");
	insist(key   != NULL, "bloom_set() is undefined for NULL keys");

	unsigned int n, i, h1, h2;
	h1 = s_djb(key, len);
	h2 = s_murmur3(key, len, h1);

	for (n = 0; n < bloom->k; n++) {
		i = (h1 + n*h2) % bloom->m;
		bloom->b[i / WORDSIZE] |= (1 << (i % WORDSIZE));
	}
}

int bloom_isset(bloom_t bloom, const char *key, size_t len)
{
	insist(bloom != NULL, "bloom_isset() is undefined for NULL bloom filters");
	insist(key   != NULL, "bloom_isset() is undefined for NULL keys");

	unsigned int n, i, h1, h2;
	h1 = s_djb(key, len);
	h2 = s_murmur3(key, len, h1);

	for (i = 1, n = 0; n < bloom->k; n++) {
		i = (h1 + n*h2) % bloom->m;
		if (!(bloom->b[i / WORDSIZE] & (1 << (i % WORDSIZE)))) {
			/* definitely not in the set */
			return 0;
		}
	}
	/* might be in the set; let's be optimistic */
	return 1;
}

int main(int argc, char **argv)
{
	const char *keys[] = {
		"A",
		"AB",
		"ABA",
		"ABBA",
		"CAR",
		"CDR",
		"CADR",
		"CADADDR",
		NULL
	};
	const char **p, **q;
	bloom_t b;

	b = bloom_new(100000, 3);
	fprintf(stderr, "false positive rate is %lf\n", bloom_fp(b));
	for (p = keys; *p; p++) {
		for (q = keys; *q; q++) {
			if (bloom_isset(b, *q, strlen(*q))) {
				printf("checking... [%s] might be in the set\n", *q);
			} else {
				printf("checking... [%s] definitely not in the set\n", *q);
			}
		}
		printf("SETTING [%s] in the filter\n", *p);
		bloom_set(b, *p, strlen(*p));
		bloom_diag(b, stderr, "");
	}

	for (q = keys; *q; q++) {
		if (bloom_isset(b, *q, strlen(*q))) {
			printf("checking... [%s] might be in the set\n", *q);
		} else {
			printf("checking... [%s] definitely not in the set\n", *q);
		}
	}
}
