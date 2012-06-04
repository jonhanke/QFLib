
// ------------------------------------------------- Modified PARI Routines ------------------------------------------------------------------


/* Minimal vectors for the integral definite quadratic form: a.
 * Result u:
 *   u[1]= Number of vectors of square norm <= BORNE
 *   u[2]= maximum norm found
 *   u[3]= list of vectors found (at most STOCKMAX)
 *
 *  If BORNE = gen_0: Minimal non-zero vectors.
 *  flag = min_ALL,   as above
 *  flag = min_FIRST, exits when first suitable vector is found.
 *  flag = min_PERF,  only compute rank of the family of v.v~ (v min.)
 *  flag = min_VECSMALL, return a t_VECSMALL of (half) the number of vectors for each norm
 *  flag = min_VECSMALL2, same but count only vectors with even norm, and shift the answer
 */


static GEN
jon_minim00(GEN a, GEN BORNE, GEN STOCKMAX, long flag)
{
  GEN x,res,p1,u,r,L,gnorme,invp,V;
  long n = lg(a), i, j, k, s, maxrank;
  pari_sp av0 = avma, av1, av, lim;
  double p,maxnorm,BOUND,*v,*y,*z,**q, eps = 0.000001;

  maxrank = 0; res = V = invp = NULL; /* gcc -Wall */


  // This was a switch...
  maxrank = itos(BORNE);
  res = vecsmall_const(maxrank, 0);
  if (flag == min_VECSMALL2) BORNE = shifti(BORNE,1);
  if (gcmp0(BORNE)) return res;
  break;


  if (n == 1)
  {

    // This was a switch...
    return res;

    res[1]=res[2]= (long)gen_0;
    res[3]=lgetg(1,t_MAT); return res;
  }

  av = avma;
  minim_alloc(n, &q, &x, &y, &z, &v);
  av1 = avma;

  u = lllgramint(a);
  if (lg(u) != n) err(talker,"not a definite form in minim00");
  a = qf_base_change(a,u,1);

  n--;
  a = mat_to_MP(a, DEFAULTPREC); r = sqred1(a);
  for (j=1; j<=n; j++)
  {
    v[j] = rtodbl(gcoeff(r,j,j));
    for (i=1; i<j; i++) q[i][j] = rtodbl(gcoeff(r,i,j));
  }

  if (gcmp0(BORNE))
  {
    double c, b = rtodbl(gcoeff(a,1,1));

    for (i=2; i<=n; i++) { c = rtodbl(gcoeff(a,i,i)); if (c < b) b = c; }
    BOUND = b+eps;
    BORNE = ground(dbltor(BOUND));
    maxnorm = -1.; /* don't update maxnorm */
  }
  else
  {
    BORNE = gfloor(BORNE);
    BOUND = gtodouble(BORNE)+eps;
    maxnorm = 0.;
  }


  s = 0; av1 = avma; lim = stack_lim(av1,1);
  k = n; y[n] = z[n] = 0;
  x[n] = (long) sqrt(BOUND/v[n]);

  for(;;x[1]--)
  {
    do
    {
      if (k>1)
      {
        long l = k-1;
	z[l] = 0;
	for (j=k; j<=n; j++) z[l] += q[l][j]*x[j];
	p = (double)x[k] + z[k];
	y[l] = y[k] + p*p*v[k];
	x[l] = (long) floor(sqrt((BOUND-y[l])/v[l])-z[l]);
        k = l;
      }
      for(;;)
      {
	p = (double)x[k] + z[k];
	if (y[k] + p*p*v[k] <= BOUND) break;
	k++; x[k]--;
      }
    }
    while (k > 1);
    if (! x[1] && y[1]<=eps) break;
    p = (double)x[1] + z[1]; p = y[1] + p*p*v[1]; /* norm(x) */
    if (maxnorm >= 0)
    {  
      if (p > maxnorm) maxnorm = p;
    }
    else
    {
      pari_sp av2 = avma;
      gnorme = ground(dbltor(p));
      if (gcmp(gnorme,BORNE) >= 0) avma = av2;
      else
      {
        BOUND=gtodouble(gnorme)+eps; s=0;
        affii(gnorme,BORNE); avma = av1;
      }
    }
    s++;

    switch(flag)
    {

      case min_VECSMALL:
	{
	  ulong norm = (ulong)(p + 0.5);
	  res[norm]++;
	}
	break;

      case min_VECSMALL2:
	{
	  ulong norm = (ulong)(p + 0.5);
	  if ((norm&1) == 0) res[norm>>1]++;
	}
	break;

    }
  }


  // This was a switch...
  avma=av; return res;


}

GEN
jon_qfrep0(GEN a, GEN borne, long flag)
{
  pari_sp av = avma;
  GEN g = jon_minim00(a, borne, gen_0, flag);
  if ((flag & 2) == 0) g = gerepileupto(av, gtovec(g));
  return g;
}






