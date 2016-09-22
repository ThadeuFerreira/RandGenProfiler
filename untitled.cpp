/* 
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and 
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and 
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/

#include <stdio.h>
//#include "mt19937-64/mt64.h"
//#include "mt19937-64/mt19937-64.c"
#include "mt19937ar.c"
#include "wellRNG512.c"
#include "wellRNG1024a.c"
#include <iostream>
#include <math.h>

using namespace std;




       


float box_muller( float m, float s)  /* normal random variate generator */
{               /* mean m, standard deviation s */
  float x1, x2, w, y1;
  static float y2;
  static int use_last = 0;

  if (use_last)           /* use value from previous call */
  {
    y1 = y2;
    use_last = 0;
  }
  else
  {
    int i = 0;
    do {
      float U1;
      float U2;
      if(i % 2){
         U1= genrand_int32()*FACTX;
         U2 = WELLRNG512a()*FACTX;
      }
      else{
        U1 = genrand_int32()*FACTX;
        U2 = WELLRNG1024a()*FACTX;
      }
      i++;
      
      x1 = 2.0 * U1 - 1.0;
      x2 = 2.0 * U2 - 1.0;
      w = x1 * x1 + x2 * x2;
      //cout << U1 << " - " << U2 << " - " << w << endl;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }

  return( m + y1 * s );
}

int main(void)
{
    int i;
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);    
    unsigned int v = 123;
    InitWELLRNG512a(&v);
    InitWELLRNG1024a(&v);
    //printf("1000 outputs of genrand64_int64()\n");
    int count = 50;
    for (i=0; i<count; i++) {
      //cout << Mersenne() << endl;
      //cout << WELLRNG512a() << endl;
      //cout << WELLRNG1024a() << endl;
      for (int j = 0; j < count; ++j)
      {
        for (int k = 0; k < count; ++k)
        {
                //box_muller(Mersenne()*FACTX,WELLRNG512a()*FACTX, 10, 2 );
               // box_muller(WELLRNG1024a()*FACTX, WELLRNG512a()*FACTX, 10, 2 );
          cout << box_muller(10,2) << endl;
                
        }
        /* code */
      }

    }


    //cout << "Finish3" << endl;
    return 0;
}


