/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */

#define WX 32
#define RX 32
#define M1X 3
#define M2X 24
#define M3X 10

#define MAT0POSX(t,v) (v^(v>>t))
#define MAT0NEGX(t,v) (v^(v<<(-(t))))
#define IdentityX(v) (v)

#define V0X            STATEX[STATEX_i                   ]
#define VM1X           STATEX[(STATEX_i+M1X) & 0x0000001fU]
#define VM2X           STATEX[(STATEX_i+M2X) & 0x0000001fU]
#define VM3X           STATEX[(STATEX_i+M3X) & 0x0000001fU]
#define VRm1X          STATEX[(STATEX_i+31) & 0x0000001fU]
#define newV0X         STATEX[(STATEX_i+31) & 0x0000001fU]
#define newV1X         STATEX[STATEX_i                   ]

#define FACTX 2.32830643653869628906e-10

static unsigned int STATEX_i = 0;
static unsigned int STATEX[RX];
static unsigned int z0X, z1X, z2X;

void InitWELLRNG1024a (unsigned int *init){
   int j;
   STATEX_i = 0;
   for (j = 0; j < RX; j++)
     STATEX[j] = init[j];
}

unsigned long WELLRNG1024a (void){
  z0X    = VRm1X;
  z1X    = IdentityX(V0X)       ^ MAT0POSX (8, VM1X);
  z2X    = MAT0NEGX (-19, VM2X) ^ MAT0NEGX(-14,VM3X);
  newV1X = z1X                 ^ z2X; 
  newV0X = MAT0NEGX (-11,z0X)   ^ MAT0NEGX(-7,z1X)    ^ MAT0NEGX(-13,z2X) ;
  STATEX_i = (STATEX_i + 31) & 0x0000001fU;
 // return ((double) STATEX[STATEX_i]  * FACTX);
  return STATEX[STATEX_i] ;
}
