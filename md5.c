#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

typedef union uwb {
    unsigned w;
    unsigned char b[4];
} WBunion;

typedef unsigned Digest[4];

unsigned f0( unsigned abcd[] ) {
    return (abcd[1] & abcd[2]) | (~abcd[1] & abcd[3]);
}

unsigned f1( unsigned abcd[] ) {
    return (abcd[3] & abcd[1]) | (~abcd[3] & abcd[2]);
}

unsigned f2( unsigned abcd[] ) {
    return abcd[1] ^ abcd[2] ^ abcd[3];
}

unsigned f3( unsigned abcd[] ) {
    return abcd[2] ^ (abcd[1] |~ abcd[3]);
}

typedef unsigned (*DgstFctn)(unsigned a[]);

unsigned *calcKs(unsigned *k)
{
    double s, pwr;
    int i;

    pwr = pow(2, 32);

    for(i = 0; i <64; i++) {
        s = fabs(sin(1+i));
        k[i] = (unsigned) (s * pwr);
    }

    return k;
}

// Rotate v left by amt bits
unsigned rol(unsigned v, short amt)
{
    unsigned msk1 = (1<<amt) -1;

    return ((v>>(32-amt)) & msk1) | ((v<<amt) & ~msk1);
}

//unsigned *md5(const char *msg, int mlen)
unsigned *md5(const char *msg)
// we CAN'T do this:
// char *md5(const char *msg)
{
    const int mlen = strlen(msg);
    static Digest h0 = {0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476};

    static DgstFctn ff[] = {&f0, &f1, &f2, &f3};

    static short M[] = {1, 5, 3, 7};
    static short O[] = {0, 1, 5, 0};
    static short rot0[] = {7, 12, 17, 22};
    static short rot1[] = {5, 9, 14, 20};
    static short rot2[] = {4, 11, 16, 23};
    static short rot3[] = {6, 10, 15, 21};
    static short *rots[] = {rot0, rot1, rot2, rot3};
    static unsigned kspace[64];
    static unsigned *k;

    static Digest h;
    Digest abcd;
    DgstFctn fctn;
    short m, o, g;
    unsigned f;
    short *rotn;

    union {
        unsigned w[16];
        char b[64];
    } mm;

    int os = 0;
    int grp, grps, q, p;
    unsigned char *msg2;

    if (k == NULL) {
        k = calcKs(kspace);
    }

    for(q=0; q<4; q++) {
        h[q] = h0[q]; // initialize
    }

    {
        grps = 1 + (mlen + 8) / 64;
        msg2 = malloc(64 * grps);
        memcpy(msg2, msg, mlen);
        msg2[mlen] = (unsigned char)0x80;
        q = mlen + 1;

        while(q < 64*grps) {
            msg2[q] = 0;
            q++;
        }

        {
            WBunion u;
            u.w = 8 * mlen;
            q -= 8;
            memcpy(msg2 + q, &u.w, 4);
        }
    }

    for(grp = 0; grp < grps; grp++) {
        memcpy(mm.b, msg2 + os, 64);

        for(q=0; q<4; q++) {
            abcd[q] = h[q];
        }

        for(p = 0; p<4; p ++) {
            fctn = ff[p];
            rotn = rots[p];
            m = M[p];
            o = O[p];

            for(q=0; q<16; q++) {
                g = (m*q + o) % 16;
                f = abcd[1] + rol(abcd[0] + fctn(abcd) + k[q + 16*p] + mm.w[g], rotn[q%4]);

                abcd[0] = abcd[3];
                abcd[3] = abcd[2];
                abcd[2] = abcd[1];
                abcd[1] = f;
            }
        }

        for(p=0; p<4; p++) {
            h[p] += abcd[p];
        }

        os += 64;
    }

    if(msg2) {
        free(msg2);
    }

    return h;
}

int main(int argc, char **argv)
{
	clock_t begin_run_time = clock();
    // defined some variable
    const int MD5_LENGTH = 32;
    const int ROOT_PROCESS = 0;
	// @TODO: CHANGE ACCORDING PASSWORD LENGTH
    const int LENGTH = 7; // max = 9
    const int MIN_RANGE = pow(10, LENGTH - 1);
    const int MAX_RANGE = pow(10, LENGTH); //max 2147483647
    const char *input = "cd7bb5142bee7b53a86ef1d4617e0599";

    MPI_Init(&argc, &argv);

    int my_rank, com_size, per_rank, c;
	char hostname[50];
	gethostname(hostname, 50);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &com_size);

    /*
    if (my_rank == ROOT_PROCESS) {
        // temporary variable
        char str_min_range[MD5_LENGTH], str_max_range[MD5_LENGTH], str_length[MD5_LENGTH];
        FILE *f = fopen("md5_digest", "rt");//rt for read lines
        // read min range
        fgets(str_min_range, MD5_LENGTH, f); // on line 1
        MIN_RANGE = atoi(str_min_range);
        printf("min: %d\n", MIN_RANGE);
        //read  max range
        fgets(str_max_range, MD5_LENGTH, f); // on line 2
        MAX_RANGE = atoi(str_max_range);
        printf("max: %d\n", MAX_RANGE);
        // read length of key
        //fgets(str_length, MD5_LENGTH, f);
        //LENGTH = atoi(str_length);
        //printf("length: %d\n", LENGTH);

        // read digest
        fgets(input, MD5_LENGTH+1, f); // on line 4
        printf("input: %s\n", input);

        fclose(f);

        // Send MIN_RANGE, MAX_RANGE, key LENGHT, md5 digest (input) to other process
        for(c=1; c<com_size;c++) {
            MPI_Send(MIN_RANGE, MD5_LENGTH, MPI_INT, c, 1, MPI_COMM_WORLD);
            MPI_Send(MAX_RANGE, MD5_LENGTH, MPI_INT, c, 2, MPI_COMM_WORLD);
            MPI_Send(input, MD5_LENGTH, MPI_CHAR, c, 3, MPI_COMM_WORLD);
        }

    } else {
        MPI_Recv(MIN_RANGE, MD5_LENGTH, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD);
        MPI_Recv(MAX_RANGE, MD5_LENGTH, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD);
        MPI_Recv(input, MD5_LENGTH, MPI_CHAR, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD);
    }
    */

    per_rank = (MAX_RANGE - MIN_RANGE) / com_size;

    int begin = MIN_RANGE + my_rank*per_rank;
    int end = MIN_RANGE + per_rank * (my_rank+1);
    printf("rank %d searching from %d to %d at %s\n", my_rank, begin, end-1, hostname);

	clock_t begin_time = clock();

    for(c= begin; c<end; c++) {
        char m[LENGTH*2];
        // @TODO: change "%0_d" with _ = LENGTH
        snprintf(m, LENGTH*2, "%07d", c);

        // --- begin encrypt
        unsigned *d = md5(m);

        WBunion u;
        char result[MD5_LENGTH];
        result[0] = '\0';
        char tmp[4];
        int i, j;

        for(i=0; i<4; i++) {
            u.w = d[i];

            for(j=0; j<4; j++) {
                snprintf(tmp, 4, "%02x", u.b[j]);
                strcat(result, tmp);
            }
        }
        //--- end encrypt

        //--- check with input
        int count = 0;
        for(i=0; i<MD5_LENGTH; i++) {
            if(input[i] != result[i]) {
                break;
            }

            count++;
        }

        if(count == MD5_LENGTH) {
			clock_t end_time = clock();
			double time_spent = (double)(end_time - begin_time) / CLOCKS_PER_SEC;

			printf("key: %s at rank %d, time: %f seconds, %d passwords scanded\n", m, my_rank, time_spent, c-begin+1);
            MPI_Abort(MPI_COMM_WORLD, MPI_SUCCESS);

            MPI_Finalize();

            return 0;
        }
    }

  	clock_t end_run_time = clock();
  	double run_time_spent = (double)(end_run_time - begin_run_time) / CLOCKS_PER_SEC;
  	printf("RUN TIME: %f seconds (please take a MAX number)\n", run_time_spent);

    printf("rank %d not found\n", my_rank);

    MPI_Finalize();

    return 0;
}
