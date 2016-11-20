/*
Kod Ÿród³owy funkcji do transformacji wspó³rzêdnych elipsoidalnych na p³askie odwzorowañ kartograficznych UTM, 1992, 2000
Krzysztof Naus, Gdynia 2014
*/
//---------------------------------------------------------------------------
#ifndef Unit_UTM_1992_2000_H
#define Unit_UTM_1992_2000_H


#include <math.h>
////////////////////////////////////////////////////////////////////////////////
//Deklaracja sta³ych
////////////////////////////////////////////////////////////////////////////////


static const double fe = 500000.0;

//Deklaracja tablicy stref rownole¿nikowych
static char cArray[] = "CDEFGHJKLMNPQRSTUVWX";

//------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
//Funkcje pomocnicze
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------

double CalculateESquared (double a, double b)
{
	return ((a * a) - (b * b)) / (a * a);
}


double CalculateE2Squared (double a, double b)
{
	return ((a * a) - (b * b)) / (b * b);
}


double denom (double es, double sphi)
{
	double sinSphi = sin (sphi);
	return sqrt (1.0 - es * (sinSphi * sinSphi));
}


double sphsr (double a, double es, double sphi)
{
	double dn = denom (es, sphi);
	return a * (1.0 - es) / (dn * dn * dn);
}


double sphsn (double a, double es, double sphi)
{
	double sinSphi = sin (sphi);
	return a / sqrt (1.0 - es * (sinSphi * sinSphi));
}


double sphtmd (double ap, double bp, double cp, double dp, double ep, double sphi)
{
	return (ap * sphi) - (bp * sin (2.0 * sphi)) + (cp * sin (4.0 * sphi))
		- (dp * sin (6.0 * sphi)) + (ep * sin (8.0 * sphi));
}

//=======================================================================
//  Funkcja s³u¿y do konwersji wspó³rzednych elipsoidalnych B, L (lat/lon) na p³askie X-northing, Y-easting odwzorowania kartograficznego UTM (dla dowolnej elipsoidy)
//=======================================================================
//      Argumenty wejœciowe i wyjœciowe:
//      --------------------------------
//       double a: d³ugoœæ du¿ej pó³osi elipsoidy odniesienia, w metrach (np. dla elipsoidy WGS 84, 6378137.0)
//       double f: sp³aszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
//       int& utmXZone: nr strefy UTM wg. podzia³u po³udnikowego (zwracane numery od 1 do 60, ka¿da strefa ma szeœæ stopni)
//       char& utmYZone: nr strefy wg. podzia³u równole¿nikowego (zwracane wartoœci: CDEFGHJKLMNPQRSTUVWX)
//       double& easting: wspó³rzêdna Y UTM, w metrach po konwersji [metry]
//       double& northing: wspó³rzêdna X UTM, w metrach po konwersji [metry]
//       double lat, double lon: wspó³rzêdne lat/lon do konwersji [stopnie]
//=======================================================================
void LatLonToUtm (double a, double f, int& utmXZone, char& utmYZone, double& easting, double& northing, double lat, double lon)
{
        // Wspó³czynnik zniekszta³cenia skali w po³udniku osiowym
        static const double ok = 0.9996;
	double recf;
	double b;
	double eSquared;
	double e2Squared;
	double tn;
	double ap;
	double bp;
	double cp;
	double dp;
	double ep;
	double olam;
	double dlam;
	double s;
	double c;
	double t;
	double eta;
	double sn;
	double tmd;
	double t1, t2, t3, t4, t5,  t6, t7, t8, t9;
	double nfn;
	if(lon <= 0.0)
          {
           utmXZone = 30 + (int)(lon / 6.0);
	  }
          else
          {
          utmXZone = 31 + (int)(lon / 6.0);
          }
        if(lat < 84.0 && lat >= 72.0)
          {
          // Specjalne zatrzymanie: strefa X ma 12 stopni od pó³nocy do po³udnia, nie 8
          utmYZone = cArray[19];
          }
          else
          {
          utmYZone = cArray[(int)((lat + 80.0) / 8.0)];
          }
	if(lat >= 84.0 || lat < -80.0)
          {
          // B³êdna wartoœæ szerokoœci geograficznej (zwracany znak gwiazdki)
           utmYZone = '*';
          }
	double latRad = lat * deg2rad;
	double lonRad = lon * deg2rad;
	recf = 1.0 / f;
	b = a * (recf - 1.0) / recf;
	eSquared = CalculateESquared (a, b);
	e2Squared = CalculateE2Squared (a, b);
	tn = (a - b) / (a + b);
	ap = a * (1.0 - tn + 5.0 * ((tn * tn) - (tn * tn * tn)) / 4.0 + 81.0 * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 64.0);
	bp = 3.0 * a * (tn - (tn * tn) + 7.0 * ((tn * tn * tn) - (tn * tn * tn * tn)) / 8.0 + 55.0 * (tn * tn * tn * tn * tn) / 64.0) / 2.0;
	cp = 15.0 * a * ((tn * tn) - (tn * tn * tn) + 3.0 * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 4.0) / 16.0;
	dp = 35.0 * a * ((tn * tn * tn) - (tn * tn * tn * tn) + 11.0 * (tn * tn * tn * tn * tn) / 16.0) / 48.0;
	ep = 315.0 * a * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 512.0;
	olam = (utmXZone * 6 - 183) * deg2rad;
	dlam = lonRad - olam;
	s = sin (latRad);
	c = cos (latRad);
	t = s / c;
	eta = e2Squared * (c * c);
	sn = sphsn (a, eSquared, latRad);
	tmd = sphtmd (ap, bp, cp, dp, ep, latRad);
	t1 = tmd * ok;
	t2 = sn * s * c * ok / 2.0;
	t3 = sn * s * (c * c * c) * ok * (5.0 - (t * t) + 9.0 * eta + 4.0 * (eta * eta)) / 24.0;
        t4 = sn * s * (c * c * c * c * c) * ok * (61.0 - 58.0 * (t * t) + (t * t * t * t) + 270.0 * eta - 330.0 * (t * t) * eta + 445.0 * (eta * eta) + 324.0 * (eta * eta * eta) - 680.0 * (t * t) * (eta * eta) + 88.0 * (eta * eta * eta * eta) - 600.0 * (t * t) * (eta * eta * eta) - 192.0 * (t * t) * (eta * eta * eta * eta)) / 720.0;
        t5 = sn * s * (c * c * c * c * c * c * c) * ok * (1385.0 - 3111.0 * (t * t) + 543.0 * (t * t * t * t) - (t * t * t * t * t * t)) / 40320.0;
	if (latRad < 0.0) nfn = 10000000.0; else nfn = 0;
	northing = nfn + t1 + (dlam * dlam) * t2 + (dlam * dlam * dlam * dlam) * t3 + (dlam * dlam * dlam * dlam * dlam * dlam) * t4 + (dlam * dlam * dlam * dlam * dlam * dlam * dlam * dlam) * t5;
	t6 = sn * c * ok;
	t7 = sn * (c * c * c) * ok * (1.0 - (t * t) + eta) / 6.0;
        t8 = sn * (c * c * c * c * c) * ok * (5.0 - 18.0 * (t * t) + (t * t * t * t) + 14.0 * eta - 58.0 * (t * t) * eta + 13.0 * (eta * eta) + 4.0 * (eta * eta * eta) - 64.0 * (t * t) * (eta * eta) - 24.0 * (t * t) * (eta * eta * eta)) / 120.0;
        t9 = sn * (c * c * c * c * c * c * c) * ok * (61.0 - 479.0 * (t * t) + 179.0 *  (t * t * t * t) - (t * t * t * t * t * t)) / 5040.0;
	easting = fe + dlam * t6 + (dlam * dlam * dlam) * t7  + (dlam * dlam * dlam * dlam * dlam) * t8 + (dlam * dlam * dlam * dlam * dlam * dlam * dlam) * t9;
	if (northing >= 9999999.0) northing = 9999999.0;
}

//=======================================================================
//  Funkcja s³u¿y do konwersji wspó³rzednych elipsoidalnych B, L (lat/lon) na p³askie X-northing, Y-easting odwzorowania kartograficznego 1992 i 2000 (dla dowolnej elipsoidy)
//=======================================================================
//      Argumenty wejœciowe i wyjœciowe:
//      --------------------------------
//       double a: d³ugoœæ du¿ej pó³osi elipsoidy odniesienia, w metrach (np. dla elipsoidy WGS 84, 6378137.0)
//       double f: sp³aszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
//       int& utmXZone: nr strefy UTM wg. podzia³u po³udnikowego (zwracane numery od 1 do 60, ka¿da strefa ma szeœæ stopni)
//       char& utmYZone: nr strefy wg. podzia³u równole¿nikowego (zwracane wartoœci: CDEFGHJKLMNPQRSTUVWX)
//       double& easting: wspó³rzêdna Y UTM, w metrach po konwersji [metry]
//       double& northing: wspó³rzêdna X UTM, w metrach po konwersji [metry]
//       double lat, double lon: wspó³rzêdne lat/lon do konwersji [stopnie]
//       int proj: odwzorowanie kartograficzne (proj = 1 odpowiada odwzorowaniu 1992, natomiast ka¿da inna odwzorowaniu 2000)
//=======================================================================
void LatLonToPUWG (double a, double f, double& easting, double& northing, double lat, double lon, int proj)
{
        // Wspó³czynnik zniekszta³cenia skli mapy w po³udniku osiowym dla odwzorowania kartograficznego 2000
        double ok =  0.999923;
        //Wspó³czynnik zniekszta³cenia skli mapy w po³udniku osiowym dla odwzorowania kartograficznego 1992
        if(proj==1) ok =  0.9993;
	double recf;
	double b;
	double eSquared;
	double e2Squared;
	double tn;
	double ap;
	double bp;
	double cp;
	double dp;
	double ep;
	double olam;
	double dlam;
	double s;
	double c;
	double t;
	double eta;
	double sn;
	double tmd;
	double t1, t2, t3, t4, t5,  t6, t7, t8, t9;
	double nfn;
        double strf;
        if(lon < 13.5 || lon > 25.5)
         {
		//B³êdna wartoœæ d³ugoœci geograficznej (zwracana wartoœæ 99999999999999)
                easting=999999999999999;
                northing=999999999999999;
                return;
	 }
         else
         {
           if(proj==1)
             {
              olam = 19.0 * deg2rad;
              strf = 0.0;
              nfn = -5300000.0;
             }
             else
             {
              nfn = 0;
              if(lon >= 13.5 && lon < 16.5)
               {
                olam = 15.0 * deg2rad;
                strf = 5000000.0;
               }

              if(lon >= 16.5 && lon < 19.5)
               {
                olam = 18.0 * deg2rad;
                strf = 6000000.0;
               }

              if(lon >= 19.5 && lon < 22.5)
               {
                olam = 21.0 * deg2rad;
                strf = 7000000.0;
               }

              if(lon >= 22.5 && lon < 25.5)
               {
                olam = 24.0 * deg2rad;
                strf = 8000000.0;
               }
             }
         }
	double latRad = lat * deg2rad;
	double lonRad = lon * deg2rad;
	recf = 1.0 / f;
	b = a * (recf - 1.0) / recf;
	eSquared = CalculateESquared (a, b);
	e2Squared = CalculateE2Squared (a, b);
	tn = (a - b) / (a + b);
	ap = a * (1.0 - tn + 5.0 * ((tn * tn) - (tn * tn * tn)) / 4.0 + 81.0 * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 64.0);
	bp = 3.0 * a * (tn - (tn * tn) + 7.0 * ((tn * tn * tn) - (tn * tn * tn * tn)) / 8.0 + 55.0 * (tn * tn * tn * tn * tn) / 64.0) / 2.0;
	cp = 15.0 * a * ((tn * tn) - (tn * tn * tn) + 3.0 * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 4.0) / 16.0;
	dp = 35.0 * a * ((tn * tn * tn) - (tn * tn * tn * tn) + 11.0 * (tn * tn * tn * tn * tn) / 16.0) / 48.0;
	ep = 315.0 * a * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 512.0;
	dlam = lonRad - olam;
	s = sin (latRad);
	c = cos (latRad);
	t = s / c;
	eta = e2Squared * (c * c);
	sn = sphsn (a, eSquared, latRad);
	tmd = sphtmd (ap, bp, cp, dp, ep, latRad);
	t1 = tmd * ok;
	t2 = sn * s * c * ok / 2.0;
	t3 = sn * s * (c * c * c) * ok * (5.0 - (t * t) + 9.0 * eta + 4.0 * (eta * eta)) / 24.0;
        t4 = sn * s * (c * c * c * c * c) * ok * (61.0 - 58.0 * (t * t) + (t * t * t * t) + 270.0 * eta - 330.0 * (t * t) * eta + 445.0 * (eta * eta) + 324.0 * (eta * eta * eta) - 680.0 * (t * t) * (eta * eta) + 88.0 * (eta * eta * eta * eta) - 600.0 * (t * t) * (eta * eta * eta) - 192.0 * (t * t) * (eta * eta * eta * eta)) / 720.0;
        t5 = sn * s * (c * c * c * c * c * c * c) * ok * (1385.0 - 3111.0 * (t * t) + 543.0 * (t * t * t * t) - (t * t * t * t * t * t)) / 40320.0;
	northing = nfn + t1 + (dlam * dlam) * t2 + (dlam * dlam * dlam * dlam) * t3 + (dlam * dlam * dlam * dlam * dlam * dlam) * t4 + (dlam * dlam * dlam * dlam * dlam * dlam * dlam * dlam) * t5;
	t6 = sn * c * ok;
	t7 = sn * (c * c * c) * ok * (1.0 - (t * t) + eta) / 6.0;
        t8 = sn * (c * c * c * c * c) * ok * (5.0 - 18.0 * (t * t) + (t * t * t * t) + 14.0 * eta - 58.0 * (t * t) * eta + 13.0 * (eta * eta) + 4.0 * (eta * eta * eta) - 64.0 * (t * t) * (eta * eta) - 24.0 * (t * t) * (eta * eta * eta)) / 120.0;
        t9 = sn * (c * c * c * c * c * c * c) * ok * (61.0 - 479.0 * (t * t) + 179.0 *  (t * t * t * t) - (t * t * t * t * t * t)) / 5040.0;
	easting = fe + strf + dlam * t6 + (dlam * dlam * dlam) * t7  + (dlam * dlam * dlam * dlam * dlam) * t8 + (dlam * dlam * dlam * dlam * dlam * dlam * dlam) * t9;// + 0.5;
}


//=======================================================================
//  Funkcja do konwersji wspó³rzêdnych p³askich X/Y UTM na elipsoidalne lat/lon (dla dowolnej elipsoidy)
//=======================================================================
//      Wymagania:
//      -------------------------------------
//      utmXZone musi byæ wartoœci¹ w garanicach od 1 do 60
//      utmYZone musi byæ jedn¹ z liter: CDEFGHJKLMNPQRSTUVWX
//      Argumenty wejœciowe i wyjœciowe:
//      ------------------------------------
//      double a: d³ugoœæ du¿ej pó³osi, w metrach (np. dla elipsoidy WGS 84, 6378137.0)
//      double f: sp³aszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
//      int utmXZone: nr strefy UTM wg. podzia³u po³udnikowego (zwracane numery od 1 do 60, ka¿da strefa ma szeœæ stopni)
//      char utmYZone: nr strefy wg. podzia³u równole¿nikowego (zwracane wartoœci: CDEFGHJKLMNPQRSTUVWX)
//      double easting, double northing: wspó³rzêdna X, Y UTM do konwersji [metry]
//      double& lat,  double& lon: wspó³rzêdne elipsoidalne lat/lon po konwersji [stopnie]
//=======================================================================
void UtmToLatLon (double a, double f, int utmXZone, char utmYZone, double easting, double northing, double& lat, double& lon)
{
        // Wspó³czynnik zniekszta³cenia skli mapy w po³udniku osiowym dla odwzorowania kartograficznego UTM
        static const double ok = 0.9996;
	double recf;
	double b;
	double eSquared;
	double e2Squared;
	double tn;
	double ap;
	double bp;
	double cp;
	double dp;
	double ep;
	double nfn;
	double tmd;
	double sr;
	double sn;
	double ftphi;
	double s;
	double c;
	double t;
	double eta;
	double de;
	double dlam;
	double olam;
	recf = 1.0 / f;
	b = a * (recf - 1) / recf;
	eSquared = CalculateESquared (a, b);
	e2Squared = CalculateE2Squared (a, b);
	tn = (a - b) / (a + b);
	ap = a * (1.0 - tn + 5.0 * ((tn * tn) - (tn * tn * tn)) / 4.0 + 81.0 * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 64.0);
	bp = 3.0 * a * (tn - (tn * tn) + 7.0 * ((tn * tn * tn) - (tn * tn * tn * tn)) / 8.0 + 55.0 * (tn * tn * tn * tn * tn) / 64.0) / 2.0;
	cp = 15.0 * a * ((tn * tn) - (tn * tn * tn) + 3.0 * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 4.0) / 16.0;
	dp = 35.0 * a * ((tn * tn * tn) - (tn * tn * tn * tn) + 11.0 * (tn * tn * tn * tn * tn) / 16.0) / 48.0;
	ep = 315.0 * a * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 512.0;
	if((utmYZone <= 'M' && utmYZone >= 'C') || (utmYZone <= 'm' && utmYZone >= 'c'))
        {
         nfn = 10000000.0;
	}
        else
        {
         nfn = 0;
	}
	tmd = (northing - nfn) / ok;
	sr = sphsr (a, eSquared, 0.0);
	ftphi = tmd / sr;
	double t10, t11, t12, t13, t14, t15, t16, t17;
	for(int i = 0; i < 5; i++)
        {
         t10 = sphtmd (ap, bp, cp, dp, ep, ftphi);
         sr = sphsr (a, eSquared, ftphi);
         ftphi = ftphi + (tmd - t10) / sr;
	}
	sr = sphsr (a, eSquared, ftphi);
	sn = sphsn (a, eSquared, ftphi);
	s = sin (ftphi);
	c = cos (ftphi);
	t = s / c;
	eta = e2Squared * (c * c);
	de = easting - fe;
	t10 = t / (2.0 * sr * sn * (ok * ok));
	t11 = t * (5.0 + 3.0 * (t * t) + eta - 4.0 * (eta * eta) - 9.0 * (t * t) * eta) / (24.0 * sr * (sn * sn * sn) * (ok * ok * ok * ok));
        t12 = t *  (61.0 + 90.0 * (t*t) + 46.0 * eta + 45.0 * (t* t * t * t) - 252.0 * (t * t) * eta - 3.0 * (eta * eta) + 100.0 * (eta * eta * eta) - 66.0 * (t * t) * (eta * eta) - 90.0 * (t * t * t * t) * eta + 88.0 * (eta * eta * eta * eta) + 225.0 * (t * t * t * t) * (eta * eta) + 84.0 * (t * t) * (eta * eta * eta) - 192.0 * (t * t) * (eta * eta * eta * eta)) / (720.0 * sr * (sn * sn * sn* sn * sn ) * (ok * ok * ok * ok * ok * ok));
        t13 = t * (1385.0 + 3633 * (t * t) + 4095.0 * (t * t * t * t) + 1575.0  * (t * t * t * t * t *t)) / (40320 * sr * (sn * sn * sn* sn * sn * sn * sn ) * (ok * ok * ok * ok * ok * ok * ok * ok));
	lat = ftphi - (de * de) * t10 + (de * de * de * de) * t11 - (de * de * de * de * de * de) * t12 + (de * de * de * de * de * de * de * de) * t13;
	t14 = 1.0 / (sn * c * ok);
	t15 = (1.0 + 2.0 * (t * t) + eta) / (6.0 * (sn * sn * sn) * c * (ok * ok * ok));
        t16 = 1.0 * (5.0 + 6.0 * eta + 28.0 * (t * t) - 3.0 * (eta * eta) + 8.0 * (t * t) * eta + 24.0 * (t * t * t * t) - 4.0 * (eta * eta * eta) + 4.0 *(t * t) * (eta * eta) + 24.0 * (t * t) * (eta * eta * eta)) / (120.0 * (sn * sn * sn * sn * sn) * c * (ok * ok * ok * ok * ok));
        t17 = 1.0 * (61.0 + 662.0 * (t * t) + 1320.0 * (t * t * t * t) + 720.0 * (t * t * t * t * t * t)) / (5040.0 * (sn * sn * sn * sn * sn * sn * sn) * c * (ok * ok * ok * ok * ok * ok * ok));
	dlam = de * t14 - (de * de * de) * t15 + (de * de * de * de * de) * t16 - (de * de * de * de * de * de * de) * t17;
	olam = (utmXZone * 6 - 183.0) * deg2rad;
	lon = olam + dlam;
	lon *= rad2deg;
	lat *= rad2deg;
}

//=======================================================================
//  Funkcja do konwersji wspó³rzêdnych p³askich X/Y odwzorowania kartograficznego 1992 i 2000 na elipsoidalne lat/lon (dla dowolnej elipsoidy)
//=======================================================================
//      double a: dl³ugoœæ du¿ej pó³si, w metrach (np. dla elipsoidy WGS 84, 6378137.0)
//      double f: sp³aszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
//      double easting, double northing: wspó³rzêdna X, Y  odwzorowania kartograficznego 1992 lub 2000 do konwersji [metry]
//      double& lat,  double& lon: wspó³rzêdne elipsoidalne lat/lon po konwersji [stopnie]
//      int proj: odwzorowanie kartograficzne (proj = 1 odpowiada odwzorowaniu 1992, natomiast ka¿da inna odwzorowaniu 2000)
//=======================================================================
void PUWGToLatLon (double a, double f, double easting, double northing, int proj, double& lat, double& lon)
{

//uklad 2000
double ok =  0.999923;
//uklad 1992 dla skal mniejszych 1:10000
if(proj==1)
 ok =  0.9993;

	double recf;
	double b;
	double eSquared;
	double e2Squared;
	double tn;
	double ap;
	double bp;
	double cp;
	double dp;
	double ep;
	double nfn;
	double tmd;
	double sr;
	double sn;
	double ftphi;
	double s;
	double c;
	double t;
	double eta;
	double de;
	double dlam;
	double olam;
        double strf;
	recf = 1.0 / f;
	b = a * (recf - 1) / recf;
	eSquared = CalculateESquared (a, b);
	e2Squared = CalculateE2Squared (a, b);
	tn = (a - b) / (a + b);
	ap = a * (1.0 - tn + 5.0 * ((tn * tn) - (tn * tn * tn)) / 4.0 + 81.0 * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 64.0);
	bp = 3.0 * a * (tn - (tn * tn) + 7.0 * ((tn * tn * tn) - (tn * tn * tn * tn)) / 8.0 + 55.0 * (tn * tn * tn * tn * tn) / 64.0) / 2.0;
	cp = 15.0 * a * ((tn * tn) - (tn * tn * tn) + 3.0 * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 4.0) / 16.0;
	dp = 35.0 * a * ((tn * tn * tn) - (tn * tn * tn * tn) + 11.0 * (tn * tn * tn * tn * tn) / 16.0) / 48.0;
	ep = 315.0 * a * ((tn * tn * tn * tn) - (tn * tn * tn * tn * tn)) / 512.0;
        if(proj==1)
          {
           olam = 19.0 * deg2rad;
           strf = 0.0;
           nfn = -5300000.0;
          }
          else
          {
           nfn = 0;
              if(easting<6000000.0 && easting>5000000.0)
               {
               strf = 5000000.0;
               olam = 15.0 * deg2rad;
               }
              if(easting<7000000.0 && easting>6000000.0)
               {
                strf = 6000000.0;
                olam = 18.0 * deg2rad;
               }
              if(easting<8000000.0 && easting>7000000.0)
               {
                strf = 7000000.0;
                olam = 21.0 * deg2rad;
               }
              if(easting<9000000.0 && easting>8000000.0)
               {
                strf = 8000000.0;
                olam = 24.0 * deg2rad;
               }
             }
	tmd = (northing - nfn) / ok;
	sr = sphsr (a, eSquared, 0.0);
	ftphi = tmd / sr;
	double t10, t11, t12, t13, t14, t15, t16, t17;
	for (int i = 0; i < 5; i++)
        {
         t10 = sphtmd (ap, bp, cp, dp, ep, ftphi);
         sr = sphsr (a, eSquared, ftphi);
         ftphi = ftphi + (tmd - t10) / sr;
	}
	sr = sphsr (a, eSquared, ftphi);
	sn = sphsn (a, eSquared, ftphi);
	s = sin (ftphi);
	c = cos (ftphi);
	t = s / c;
	eta = e2Squared * (c * c);
	de = easting - fe - strf;
	t10 = t / (2.0 * sr * sn * (ok * ok));
	t11 = t * (5.0 + 3.0 * (t * t) + eta - 4.0 * (eta * eta) - 9.0 * (t * t) * eta) / (24.0 * sr * (sn * sn * sn) * (ok * ok * ok * ok));
        t12 = t *  (61.0 + 90.0 * (t*t) + 46.0 * eta + 45.0 * (t* t * t * t) - 252.0 * (t * t) * eta - 3.0 * (eta * eta) + 100.0 * (eta * eta * eta) - 66.0 * (t * t) * (eta * eta) - 90.0 * (t * t * t * t) * eta + 88.0 * (eta * eta * eta * eta) + 225.0 * (t * t * t * t) * (eta * eta) + 84.0 * (t * t) * (eta * eta * eta) - 192.0 * (t * t) * (eta * eta * eta * eta)) / (720.0 * sr * (sn * sn * sn* sn * sn ) * (ok * ok * ok * ok * ok * ok));
        t13 = t * (1385.0 + 3633 * (t * t) + 4095.0 * (t * t * t * t) + 1575.0  * (t * t * t * t * t *t)) / (40320 * sr * (sn * sn * sn* sn * sn * sn * sn ) * (ok * ok * ok * ok * ok * ok * ok * ok));
	lat = ftphi - (de * de) * t10 + (de * de * de * de) * t11 - (de * de * de * de * de * de) * t12 + (de * de * de * de * de * de * de * de) * t13;
	t14 = 1.0 / (sn * c * ok);
	t15 = (1.0 + 2.0 * (t * t) + eta) / (6.0 * (sn * sn * sn) * c * (ok * ok * ok));
        t16 = 1.0 * (5.0 + 6.0 * eta + 28.0 * (t * t) - 3.0 * (eta * eta) + 8.0 * (t * t) * eta + 24.0 * (t * t * t * t) - 4.0 * (eta * eta * eta) + 4.0 *(t * t) * (eta * eta) + 24.0 * (t * t) * (eta * eta * eta)) / (120.0 * (sn * sn * sn * sn * sn) * c * (ok * ok * ok * ok * ok));
        t17 = 1.0 * (61.0 + 662.0 * (t * t) + 1320.0 * (t * t * t * t) + 720.0 * (t * t * t * t * t * t)) / (5040.0 * (sn * sn * sn * sn * sn * sn * sn) * c * (ok * ok * ok * ok * ok * ok * ok));
	dlam = de * t14 - (de * de * de) * t15 + (de * de * de * de * de) * t16 - (de * de * de * de * de * de * de) * t17;
	lon = olam + dlam;
	lon *= rad2deg;
	lat *= rad2deg;
}







//------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
// Funkcje zasadnicze
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------

//==============================================================================
// Funkcja do konwersji wspó³rzêdnych X/Y 1992 lub 2000 na lat/lon elipsoide WGS 84
//==============================================================================
void PUWGToLatLonWGS84(double easting, double northing, int proj, double& lat, double& lon)
    {
    PUWGToLatLon (6378137.0, 1 / 298.257223563, easting, northing, proj, lat, lon);
    }

//==============================================================================
// Funkcja do konwersji wspó³rzêdnych lat/lon WGS 84 na X/Y UTM
//==============================================================================
void LatLonToUtmWGS84(int& utmXZone, char& utmYZone, double& easting, double& northing, double lat, double lon)
   {
   LatLonToUtm (6378137.0, 1 / 298.257223563, utmXZone, utmYZone, easting, northing, lat, lon);
   }

//==============================================================================
// Funkcja do konwersji wspó³rzêdnych lat/lon WGS 84 na X/Y 1992 lub 2000
//==============================================================================
void LatLonToPUWGWGS84(double& easting, double& northing, double lat, double lon, int proj)
   {
   //proj = 1 - dla odwzorowania kartograficznego 1992, ka¿da inna wartoœæ dla 2000
   LatLonToPUWG (6378137.0, 1 / 298.257223563, easting, northing, lat, lon, proj);
   }

//==============================================================================
//  Funkcja do konwersji wspó³rzêdnych X/Y UTM na lat/lon elipsoidalne WGS 84
//==============================================================================
void UtmToLatLonWGS84(int utmXZone, char utmYZone, double easting, double northing, double& lat, double& lon)
   {
   UtmToLatLon (6378137.0, 1 / 298.257223563, utmXZone, utmYZone, easting, northing, lat, lon);
   }

//==============================================================================
//  Funkcja do konwersji wspó³rzêdnych X/Y 1992 lub 2000 na lat/lon elipsoidalne WGS 84
//==============================================================================
void PUWGToLatLonWGS84 (double easting, double northing, int proj, double& lat, double& lon)
   {
   PUWGToLatLon(6378137.0, 1 / 298.257223563, easting, northing, proj, lat, lon);
   }

#endif
