<?php

/**
 * Kod źródłowy funkcji do transformacji wspórzędnych elipsoidalnych na płaskie
 * odwzorowań kartograficznych UTM, 1992, 2000
 * @author Krzysztof Naus, Gdynia 2014
 * @source http://ecdis.republika.pl/
 * @author Piotr Knapik (porting to PHP), 2016
 */

define('DEG2RAD', 0.01745329);
define('RAD2DEG', 57.29577951);

class UtmConverter
{

    private $fe = 500000.0;

    // Deklaracja tablicy stref równoleżnikowych
    private $cArray = "CDEFGHJKLMNPQRSTUVWX";

    /**
     * @param double $a
     * @param double $b
     * @return double
     */
    private function CalculateESquared($a, $b)
    {
        return (($a * $a) - ($b * $b)) / ($a * $a);
    }

    /**
     * @param double $a
     * @param double $b
     * @return double
     */
    private function CalculateE2Squared($a, $b)
    {
        return (($a * $a) - ($b * $b)) / ($b * $b);
    }

    /**
     * @param double $es
     * @param double $sphi
     * @return double
     */
    private function denom($es, $sphi)
    {
        /** @var double $sinSphi */
        $sinSphi = sin($sphi);
        return sqrt(1.0 - $es * ($sinSphi * $sinSphi));
    }

    /**
     * @param double $a
     * @param double $es
     * @param double $sphi
     * @return double
     */
    private function sphsr($a, $es, $sphi)
    {
        /** @var double $dn */
        $dn = $this->denom($es, $sphi);
        return $a * (1.0 - $es) / ($dn * $dn * $dn);
    }

    /**
     * @param double $a
     * @param double $es
     * @param double $sphi
     * @return double
     */
    private function sphsn($a, $es, $sphi)
    {
        /** @var double $sinSphi */
        $sinSphi = sin($sphi);
        return $a / sqrt(1.0 - $es * ($sinSphi * $sinSphi));
    }

    /**
     * @param double $ap
     * @param double $bp
     * @param double $cp
     * @param double $dp
     * @param double $ep
     * @param double $sphi
     * @return double
     */
    private function sphtmd($ap, $bp, $cp, $dp, $ep, $sphi)
    {
        return ($ap * $sphi) - ($bp * sin(2.0 * $sphi)) + ($cp * sin(4.0 * $sphi)) - ($dp * sin(6.0 * $sphi)) + ($ep * sin(8.0 * $sphi));
    }

    /**
     * Funkcja konwertująca współrzędne elipsoidalne B, L (lat/lon)
     * na płaskie X-northing, Y-easting odwzorowania kartograficznego UTM (dla dowolnej elipsoidy)
     *
     * @param double $a Długość dużej półosi elipsoidy odniesienia, w metrach (np. dla elipsoidy WGS 84, 6378137.0)
     * @param double $f Spłaszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
     * @param int $utmXZone nr strefy UTM wg. podziału południkowego (zwracane numery od 1 do 60, każda strefa ma sześć stopni)
     * @param string $utmYZone nr strefy wg. podziału równoleżnikowego (zwracane wartości: CDEFGHJKLMNPQRSTUVWX)
     * @param double $easting wspłrzędna Y UTM, w metrach po konwersji [metry]
     * @param double $northing wspłrzędna X UTM, w metrach po konwersji [metry]
     * @param double $lat wspórzędne lat/lon do konwersji [stopnie]
     * @param double $lon wspórzędne lat/lon do konwersji [stopnie]
     * @void
     */
    public function LatLonToUtm($a, $f, &$utmXZone, &$utmYZone, &$easting, &$northing, $lat, $lon)
    {
        /** @var double $ok Wspóczynnik zniekształcenia skali w południku osiowym */
        $ok = 0.9996;

        if ($lon <= 0.0) {
            $utmXZone = 30 + (int)($lon / 6.0);
        } else {
            $utmXZone = 31 + (int)($lon / 6.0);
        }

        if ($lat < 84.0 && $lat >= 72.0) {
            // Specjalne zatrzymanie: strefa X ma 12 stopni od północy do południa, nie 8
            $utmYZone = $this->cArray[19];
        } else {
            $utmYZone = $this->cArray[(int)(($lat + 80.0) / 8.0)];
        }
        if ($lat >= 84.0 || $lat < -80.0) {
            // Błędna wartość szerokości geograficznej (zwracany znak gwiazdki)
            $utmYZone = '*';
        }

        /** @var double $latRad */
        $latRad = $lat * DEG2RAD;
        /** @var double $lonRad */
        $lonRad = $lon * DEG2RAD;
        /** @var double $recf */
        $recf = 1.0 / $f;
        /** @var double $b */
        $b = $a * ($recf - 1.0) / $recf;
        /** @var double $eSquared */
        $eSquared = $this->CalculateESquared($a, $b);
        /** @var double $e2Squared */
        $e2Squared = $this->CalculateE2Squared($a, $b);
        /** @var double $tn */
        $tn = ($a - $b) / ($a + $b);
        /** @var double $ap */
        $ap = $a * (1.0 - $tn + 5.0 * (($tn * $tn) - ($tn * $tn * $tn)) / 4.0 + 81.0 * (($tn * $tn * $tn * $tn) - ($tn * $tn * $tn * $tn * $tn)) / 64.0);
        /** @var double $bp */
        $bp = 3.0 * $a * ($tn - ($tn * $tn) + 7.0 * (($tn * $tn * $tn) - ($tn * $tn * $tn * $tn)) / 8.0 + 55.0 * ($tn * $tn * $tn * $tn * $tn) / 64.0) / 2.0;
        /** @var double $cp */
        $cp = 15.0 * $a * (($tn * $tn) - ($tn * $tn * $tn) + 3.0 * (($tn * $tn * $tn * $tn) - ($tn * $tn * $tn * $tn * $tn)) / 4.0) / 16.0;
        /** @var double $dp */
        $dp = 35.0 * $a * (($tn * $tn * $tn) - ($tn * $tn * $tn * $tn) + 11.0 * ($tn * $tn * $tn * $tn * $tn) / 16.0) / 48.0;
        /** @var double $ep */
        $ep = 315.0 * $a * (($tn * $tn * $tn * $tn) - ($tn * $tn * $tn * $tn * $tn)) / 512.0;
        /** @var double $olam */
        $olam = ($utmXZone * 6 - 183) * DEG2RAD;
        /** @var double $dlam */
        $dlam = $lonRad - $olam;
        /** @var double $s */
        $s = sin($latRad);
        /** @var double $c */
        $c = cos($latRad);
        /** @var double $t */
        $t = $s / $c;
        /** @var double $eta */
        $eta = $e2Squared * ($c * $c);
        /** @var double $sn */
        $sn = $this->sphsn($a, $eSquared, $latRad);
        /** @var double $tmd */
        $tmd = $this->sphtmd($ap, $bp, $cp, $dp, $ep, $latRad);
        /** @var double $t1 */
        $t1 = $tmd * $ok;
        /** @var double $t2 */
        $t2 = $sn * $s * $c * $ok / 2.0;
        /** @var double $t3 */
        $t3 = $sn * $s * pow($c, 3) * $ok * (5.0 - pow($t, 2) + 9.0 * $eta + 4.0 * pow($eta, 2)) / 24.0;
        /** @var double $t4 */
        $t4 = $sn * $s * pow($c, 5) * $ok * (
                61.0
                - 58.0 * pow($t, 2) + pow($t, 4)
                + 270.0 * $eta
                - 330.0 * pow($t, 2) * $eta
                + 445.0 * pow($eta, 2)
                + 324.0 * pow($eta, 3)
                - 680.0 * pow($t, 2) * pow($eta, 2)
                + 88.0 * pow($eta, 4)
                - 600.0 * pow($t, 2) * pow($eta, 3)
                - 192.0 * pow($t, 2) * pow($eta, 4)
            ) / 720.0;
        /** @var double $t5 */
        $t5 = $sn * $s * pow($c, 7) * $ok * (1385.0 - 3111.0 * pow($t, 2) + 543.0 * pow($t, 4) - pow($t, 6)) / 40320.0;
        if ($latRad < 0.0) {
            /** @var double $nfn */
            $nfn = 10000000.0;
        } else {
            $nfn = 0;
        }
        $northing = $nfn
            + $t1
            + pow($dlam, 2) * $t2
            + pow($dlam, 4) * $t3
            + pow($dlam, 6) * $t4
            + pow($dlam, 8) * $t5;
        /** @var double $t6 */
        $t6 = $sn * $c * $ok;
        /** @var double $t7 */
        $t7 = $sn * pow($c, 3) * $ok * (1.0 - pow($t, 2) + $eta) / 6.0;
        /** @var double $t8 */
        $t8 = $sn * pow($c, 5) * $ok * (
                5.0
                - 18.0 * pow($t, 2)
                + pow($t, 4) + 14.0 * $eta
                - 58.0 * pow($t, 2) * $eta
                + 13.0 * pow($eta, 2)
                + 4.0 * pow($eta, 3)
                - 64.0 * pow($t, 2) * pow($eta, 2)
                - 24.0 * pow($t, 2) * pow($eta, 3)
            ) / 120.0;
        /** @var double $t9 */
        $t9 = $sn * pow($c, 7) * $ok * (61.0 - 479.0 * pow($t, 2) + 179.0 * pow($t, 4) - pow($t, 6)) / 5040.0;
        $easting = $this->fe + $dlam * $t6 + pow($dlam, 3) * $t7 + pow($dlam, 5) * $t8 + pow($dlam, 7) * $t9;
        if ($northing >= 9999999.0) {
            $northing = 9999999.0;
        }
    }

    /**
     * Funkcja konwertująca współrzędne elipsoidalne B, L (lat/lon)
     * na płaskie X-northing, Y-easting odwzorowania kartograficznego 1992 i 2000 (dla dowolnej elipsoidy)
     *
     * @param double $a Długość dużej półosi elipsoidy odniesienia, w metrach (np. dla elipsoidy WGS 84, 6378137.0)
     * @param double $f Spłaszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
     * @param double $easting wspłrzędna Y UTM, w metrach po konwersji [metry]
     * @param double $northing wspłrzędna X UTM, w metrach po konwersji [metry]
     * @param double $lat wspórzędne lat/lon do konwersji [stopnie]
     * @param double $lon wspórzędne lat/lon do konwersji [stopnie]
     * @param int $proj Odwzorowanie kartograficzne (proj = 1 odpowiada odwzorowaniu 1992, natomiast każda inna odwzorowaniu 2000)
     * @void
     */
    public function LatLonToPUWG($a, $f, &$easting, &$northing, $lat, $lon, $proj)
    {
        /** @var double $ok Wspóczynnik zniekształcenia skli mapy w po�udniku osiowym dla odwzorowania kartograficznego 2000 */
        $ok = 0.999923;
        if ($proj == 1) {
            // Wspóczynnik zniekształcenia skli mapy w południku osiowym dla odwzorowania kartograficznego 1992
            $ok = 0.9993;
        }

        if ($lon < 13.5 || $lon > 25.5) {
            //Błędna wartość długości geograficznej (zwracana wartość 99999999999999)
            $easting = 999999999999999;
            $northing = 999999999999999;
            return;
        }

        // defaults for $proj == 1
        /** @var double $olam */
        $olam = 19.0 * DEG2RAD;
        /** @var double $strf */
        $strf = 0.0;
        /** @var double $nfn */
        $nfn = -5300000.0;

        if ($proj != 1) {
            $nfn = 0;
            if ($lon >= 13.5 && $lon < 16.5) {
                $olam = 15.0 * DEG2RAD;
                $strf = 5000000.0;
            }

            if ($lon >= 16.5 && $lon < 19.5) {
                $olam = 18.0 * DEG2RAD;
                $strf = 6000000.0;
            }

            if ($lon >= 19.5 && $lon < 22.5) {
                $olam = 21.0 * DEG2RAD;
                $strf = 7000000.0;
            }

            if ($lon >= 22.5 && $lon < 25.5) {
                $olam = 24.0 * DEG2RAD;
                $strf = 8000000.0;
            }
        }

        /** @var double $latRad */
        $latRad = $lat * DEG2RAD;
        /** @var double $lonRad */
        $lonRad = $lon * DEG2RAD;
        /** @var double $recf */
        $recf = 1.0 / $f;
        /** @var double $b */
        $b = $a * ($recf - 1.0) / $recf;
        /** @var double $eSquared */
        $eSquared = $this->CalculateESquared($a, $b);
        /** @var double $e2Squared */
        $e2Squared = $this->CalculateE2Squared($a, $b);
        /** @var double $tn */
        $tn = ($a - $b) / ($a + $b);
        /** @var double $ap */
        $ap = $a * (1.0 - $tn + 5.0 * (pow($tn, 2) - pow($tn, 3)) / 4.0 + 81.0 * (pow($tn, 4) - pow($tn, 5)) / 64.0);
        /** @var double $bp */
        $bp = 3.0 * $a * ($tn - pow($tn, 2) + 7.0 * (pow($tn, 3) - pow($tn, 4)) / 8.0 + 55.0 * pow($tn, 5) / 64.0) / 2.0;
        /** @var double $cp */
        $cp = 15.0 * $a * (pow($tn, 2) - pow($tn, 3) + 3.0 * (pow($tn, 4) - pow($tn, 5)) / 4.0) / 16.0;
        /** @var double $dp */
        $dp = 35.0 * $a * (pow($tn, 3) - pow($tn, 4) + 11.0 * pow($tn, 5) / 16.0) / 48.0;
        /** @var double $ep */
        $ep = 315.0 * $a * (pow($tn, 4) - pow($tn, 5)) / 512.0;
        /** @var double $dlam */
        $dlam = $lonRad - $olam;
        /** @var double $s */
        $s = sin($latRad);
        /** @var double $c */
        $c = cos($latRad);
        /** @var double $t */
        $t = $s / $c;
        /** @var double $eta */
        $eta = $e2Squared * ($c * $c);
        /** @var double $sn */
        $sn = $this->sphsn($a, $eSquared, $latRad);
        /** @var double $tmd */
        $tmd = $this->sphtmd($ap, $bp, $cp, $dp, $ep, $latRad);
        /** @var double $t1 */
        $t1 = $tmd * $ok;
        /** @var double $t2 */
        $t2 = $sn * $s * $c * $ok / 2.0;
        /** @var double $t3 */
        $t3 = $sn * $s * pow($c, 3) * $ok * (5.0 - pow($t, 2) + 9.0 * $eta + 4.0 * pow($eta, 2)) / 24.0;
        /** @var double $t4 */
        $t4 = $sn * $s * pow($c, 5) * $ok * (61.0 - 58.0 * pow($t, 2)
                + pow($t, 4) + 270.0 * $eta - 330.0 * pow($t, 2) * $eta
                + 445.0 * pow($eta, 2)
                + 324.0 * pow($eta, 3)
                - 680.0 * pow($t, 2) * pow($eta, 2)
                + 88.0 * pow($eta, 4)
                - 600.0 * pow($t, 2) * pow($eta, 3)
                - 192.0 * pow($t, 2) * pow($eta, 4)) / 720.0;
        /** @var double $t5 */
        $t5 = $sn * $s * pow($c, 7) * $ok * (1385.0 - 3111.0 * pow($t, 2) + 543.0 * pow($t, 4) - pow($t, 6)) / 40320.0;
        $northing = $nfn + $t1 + pow($dlam, 2) * $t2 + pow($dlam, 4) * $t3 + pow($dlam, 6) * $t4 + pow($dlam, 8) * $t5;
        /** @var double $t6 */
        $t6 = $sn * $c * $ok;
        /** @var double $t7 */
        $t7 = $sn * pow($c, 3) * $ok * (1.0 - pow($t, 2) + $eta) / 6.0;
        /** @var double $t8 */
        $t8 = $sn * pow($c, 5) * $ok * (5.0 - 18.0 * pow($t, 2)
                + pow($t, 4) + 14.0 * $eta - 58.0 * pow($t, 2) * $eta
                + 13.0 * pow($eta, 2) + 4.0 * pow($eta, 3)
                - 64.0 * pow($t, 2) * pow($eta, 2)
                - 24.0 * pow($t, 2) * pow($eta, 3)) / 120.0;
        /** @var double $t9 */
        $t9 = $sn * pow($c, 7) * $ok * (61.0 - 479.0 * pow($t, 2) + 179.0 * pow($t, 4) - pow($t, 6)) / 5040.0;
        $easting = $this->fe + $strf + $dlam * $t6 + pow($dlam, 3) * $t7 + pow($dlam, 5) * $t8 + pow($dlam, 7) * $t9;// + 0.5;
    }

    /**
     * Funkcja do konwersji wspórzłdnych płaskich X/Y UTM na elipsoidalne lat/lon (dla dowolnej elipsoidy)
     *
     * @param double $a długość dużej półosi, w metrach (np. dla elipsoidy WGS 84, 6378137.0)
     * @param double $f spłaszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
     * @param int $utmXZone nr strefy UTM wg. podziału południkowego (zwracane numery od 1 do 60, ka�da strefa ma sze�� stopni)
     * @param string $utmYZone nr strefy wg. podziału równoleżnikowego (zwracane wartości: CDEFGHJKLMNPQRSTUVWX)
     * @param double $easting wspłrzędna X UTM do konwersji [metry]
     * @param double $northing wspłrzędna Y UTM do konwersji [metry]
     * @param double &$lat współrzędna elipsoidalna lat po konwersji [stopnie]
     * @param double &$lon współrzędna elipsoidalna lon po konwersji [stopnie]
     * @void
     */
    private function UtmToLatLon($a, $f, $utmXZone, $utmYZone, $easting, $northing, &$lat, &$lon)
    {
        /** @var double $ok Wspóczynnik zniekształcenia skli mapy w południku osiowym dla odwzorowania kartograficznego UTM */
        $ok = 0.9996;
        /** @var double $recf */
        $recf = 1.0 / $f;
        /** @var double $b */
        $b = $a * ($recf - 1) / $recf;
        /** @var double $eSquared */
        $eSquared = $this->CalculateESquared($a, $b);
        /** @var double $e2Squared */
        $e2Squared = $this->CalculateE2Squared($a, $b);
        /** @var double $tn */
        $tn = ($a - $b) / ($a + $b);
        /** @var double $ap */
        $ap = $a * (1.0 - $tn + 5.0 * (pow($tn, 2) - pow($tn, 3)) / 4.0 + 81.0 * (pow($tn, 4) - pow($tn, 5)) / 64.0);
        /** @var double $bp */
        $bp = 3.0 * $a * ($tn - pow($tn, 2) + 7.0 * (pow($tn, 3) - pow($tn, 4)) / 8.0 + 55.0 * pow($tn, 5) / 64.0) / 2.0;
        /** @var double $cp */
        $cp = 15.0 * $a * (pow($tn, 2) - pow($tn, 3) + 3.0 * (pow($tn, 4) - pow($tn, 5)) / 4.0) / 16.0;
        /** @var double $dp */
        $dp = 35.0 * $a * (pow($tn, 3) - pow($tn, 4) + 11.0 * pow($tn, 5) / 16.0) / 48.0;
        /** @var double $ep */
        $ep = 315.0 * $a * (pow($tn, 4) - pow($tn, 5)) / 512.0;

        if (($utmYZone <= 'M' && $utmYZone >= 'C') || ($utmYZone <= 'm' && $utmYZone >= 'c')) {
            $nfn = 10000000.0;
        } else {
            $nfn = 0;
        }

        /** @var double tmd */
        $tmd = ($northing - $nfn) / $ok;
        /** @var double sr */
        $sr = $this->sphsr($a, $eSquared, 0.0);
        /** @var double $ftphi */
        $ftphi = $tmd / $sr;
        for ($i = 0; $i < 5; $i++) {
            $t10 = $this->sphtmd($ap, $bp, $cp, $dp, $ep, $ftphi);
            $sr = $this->sphsr($a, $eSquared, $ftphi);
            $ftphi = $ftphi + ($tmd - $t10) / $sr;
        }
        $sr = $this->sphsr($a, $eSquared, $ftphi);
        /** @var double $sn */
        $sn = $this->sphsn($a, $eSquared, $ftphi);
        /** @var double $s */
        $s = sin($ftphi);
        /** @var double $c */
        $c = cos($ftphi);
        /** @var double $t */
        $t = $s / $c;
        /** @var double $eta */
        $eta = $e2Squared * pow($c, 2);
        $de = $easting - $this->fe;
        /** @var double $t10 */
        $t10 = $t / (2.0 * $sr * $sn * pow($ok, 2));
        /** @var double $t11 */
        $t11 = $t * (5.0 + 3.0 * pow($t, 2) + $eta - 4.0 * pow($eta, 2) - 9.0 * pow($t, 2) * $eta) / (24.0 * $sr * pow($sn, 3) * pow($ok, 4));
        /** @var double $t12 */
        $t12 = $t * (
                61.0
                + 90.0 * pow($t, 2)
                + 46.0 * $eta
                + 45.0 * pow($t, 4)
                - 252.0 * pow($t, 2) * $eta
                - 3.0 * pow($eta, 2)
                + 100.0 * pow($eta, 3)
                - 66.0 * pow($t, 2) * pow($eta, 2)
                - 90.0 * pow($t, 4) * $eta
                + 88.0 * pow($eta, 4)
                + 225.0 * pow($t, 4) * pow($eta, 2)
                + 84.0 * pow($t, 2) * pow($eta, 3)
                - 192.0 * pow($t, 2) * pow($eta, 4)
            ) / (720.0 * $sr * pow($sn, 5) * pow($ok, 6));
        /** @var double $t13 */
        $t13 = $t * (1385.0 + 3633 * pow($t, 2) + 4095.0 * pow($t, 4) + 1575.0 * pow($t, 6)) / (40320 * $sr * pow($sn, 7) * pow($ok, 8));
        $lat = $ftphi - pow($de, 2) * $t10 + pow($de, 4) * $t11 - pow($de, 6) * $t12 + pow($de, 8) * $t13;
        /** @var double $t14 */
        $t14 = 1.0 / ($sn * $c * $ok);
        /** @var double $t15 */
        $t15 = (1.0 + 2.0 * pow($t, 2) + $eta) / (6.0 * pow($sn, 3) * $c * pow($ok, 3));
        /** @var double $t16 */
        $t16 = 1.0 * (
                5.0
                + 6.0 * $eta
                + 28.0 * pow($t, 2)
                - 3.0 * pow($eta, 2)
                + 8.0 * pow($t, 2) * $eta
                + 24.0 * pow($t, 4)
                - 4.0 * pow($eta, 3)
                + 4.0 * pow($t, 2) * pow($eta, 2)
                + 24.0 * pow($t, 2) * pow($eta, 3)
            ) / (120.0 * pow($sn, 5) * $c * pow($ok, 5));
        /** @var double $t17 */
        $t17 = 1.0 * (61.0 + 662.0 * pow($t, 2) + 1320.0 * pow($t, 4) + 720.0 * pow($t, 6)) / (5040.0 * pow($sn, 7) * $c * pow($ok, 7));
        /** @var double $dlam */
        $dlam = $de * $t14 - pow($de, 3) * $t15 + pow($de, 5) * $t16 - pow($de, 7) * $t17;
        /** @var double $olam */
        $olam = ($utmXZone * 6 - 183.0) * DEG2RAD;
        $lon = $olam + $dlam;
        $lon *= RAD2DEG;
        $lat *= RAD2DEG;
    }

    /**
     * Funkcja do konwersji współrzędnych płaskich X/Y odwzorowania kartograficznego 1992 i 2000 na elipsoidalne lat/lon (dla dowolnej elipsoidy)
     *
     * @param double $a długość dużej półsi, w metrach (np. dla elipsoidy WGS 84, 6378137.0)
     * @param double $f spłaszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
     * @param double $easting spłaszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
     * @param double $northing spłaszczenie elipsoidalne (np. dla elipsoidy WGS 84, 1 / 298.257223563)
     * @param int $proj odwzorowanie kartograficzne (proj = 1 odpowiada odwzorowaniu 1992, natomiast ka�da inna odwzorowaniu 2000)
     * @param double &$lat współrzędna elipsoidalna lat po konwersji [stopnie]
     * @param double &$lon współrzędna elipsoidalna lon po konwersji [stopnie]
     */
    private function PUWGToLatLon($a, $f, $easting, $northing, $proj, &$lat, &$lon)
    {
        //uklad 2000
        /** @var double $ok */
        $ok = 0.999923;
        if ($proj == 1) {
            //uklad 1992 dla skal mniejszych 1:10000
            $ok = 0.9993;
        }

        /** @var double $recf */
        $recf = 1.0 / $f;
        /** @var double $b */
        $b = $a * ($recf - 1) / $recf;
        /** @var double $eSquared */
        $eSquared = $this->CalculateESquared($a, $b);
        /** @var double $e2Squared */
        $e2Squared = $this->CalculateE2Squared($a, $b);
        /** @var double $tn */
        $tn = ($a - $b) / ($a + $b);
        /** @var double $ap */
        $ap = $a * (1.0 - $tn + 5.0 * (pow($tn, 2) - pow($tn, 3)) / 4.0 + 81.0 * (pow($tn, 4) - pow($tn, 5)) / 64.0);
        /** @var double $bp */
        $bp = 3.0 * $a * ($tn - pow($tn, 2) + 7.0 * (pow($tn, 3) - pow($tn, 4)) / 8.0 + 55.0 * pow($tn, 5) / 64.0) / 2.0;
        /** @var double $cp */
        $cp = 15.0 * $a * (pow($tn, 2) - pow($tn, 3) + 3.0 * (pow($tn, 4) - pow($tn, 5)) / 4.0) / 16.0;
        /** @var double $dp */
        $dp = 35.0 * $a * (pow($tn, 3) - pow($tn, 4) + 11.0 * pow($tn, 5) / 16.0) / 48.0;
        /** @var double $ep */
        $ep = 315.0 * $a * (pow($tn, 4) - pow($tn, 5)) / 512.0;

        // defaults for $proj == 1
        /** @var double $olam */
        $olam = 19.0 * DEG2RAD;
        /** @var double $strf */
        $strf = 0.0;
        /** @var double $nfn */
        $nfn = -5300000.0;

        if ($proj != 1) {
            $nfn = 0;
            if ($easting < 6000000.0 && $easting > 5000000.0) {
                $strf = 5000000.0;
                $olam = 15.0 * DEG2RAD;
            }
            if ($easting < 7000000.0 && $easting > 6000000.0) {
                $strf = 6000000.0;
                $olam = 18.0 * DEG2RAD;
            }
            if ($easting < 8000000.0 && $easting > 7000000.0) {
                $strf = 7000000.0;
                $olam = 21.0 * DEG2RAD;
            }
            if ($easting < 9000000.0 && $easting > 8000000.0) {
                $strf = 8000000.0;
                $olam = 24.0 * DEG2RAD;
            }
        }
        /** @var double $tmd */
        $tmd = ($northing - $nfn) / $ok;
        /** @var double $sr */
        $sr = $this->sphsr($a, $eSquared, 0.0);
        /** @var double $ftphi */
        $ftphi = $tmd / $sr;
        for ($i = 0; $i < 5; $i++) {
            $t10 = $this->sphtmd($ap, $bp, $cp, $dp, $ep, $ftphi);
            $sr = $this->sphsr($a, $eSquared, $ftphi);
            $ftphi = $ftphi + ($tmd - $t10) / $sr;
        }
        $sr = $this->sphsr($a, $eSquared, $ftphi);
        /** @var double $sn */
        $sn = $this->sphsn($a, $eSquared, $ftphi);
        /** @var double $s */
        $s = sin($ftphi);
        /** @var double $c */
        $c = cos($ftphi);
        /** @var double $t */
        $t = $s / $c;
        /** @var double $eta */
        $eta = $e2Squared * pow($c, 2);
        /** @var double $de */
        $de = $easting - $this->fe - $strf;
        /** @var double $t10 */
        $t10 = $t / (2.0 * $sr * $sn * pow($ok, 2));
        /** @var double $t11 */
        $t11 = $t * (5.0 + 3.0 * pow($t, 2) + $eta - 4.0 * pow($eta, 2) - 9.0 * pow($t, 2) * $eta) / (24.0 * $sr * pow($sn, 3) * pow($ok, 4));
        /** @var double $t12 */
        $t12 = $t * (
                61.0
                + 90.0 * pow($t, 2)
                + 46.0 * $eta
                + 45.0 * pow($t, 4)
                - 252.0 * pow($t, 2) * $eta
                - 3.0 * pow($eta, 2)
                + 100.0 * pow($eta, 3)
                - 66.0 * pow($t, 2) * pow($eta, 2)
                - 90.0 * pow($t, 4) * $eta
                + 88.0 * pow($eta, 4)
                + 225.0 * pow($t, 4) * pow($eta, 2)
                + 84.0 * pow($t, 2) * pow($eta, 3)
                - 192.0 * pow($t, 2) * pow($eta, 4)
            ) / (720.0 * $sr * pow($sn, 5) * pow($ok, 6));
        /** @var double $t13 */
        $t13 = $t * (1385.0 + 3633 * pow($t, 2) + 4095.0 * pow($t, 4) + 1575.0 * pow($t, 6)) / (40320 * $sr * pow($sn, 7) * pow($ok, 8));
        $lat = $ftphi - pow($de, 2) * $t10 + pow($de, 4) * $t11 - pow($de, 6) * $t12 + pow($de, 8) * $t13;
        /** @var double $t14 */
        $t14 = 1.0 / ($sn * $c * $ok);
        /** @var double $t15 */
        $t15 = (1.0 + 2.0 * pow($t, 2) + $eta) / (6.0 * pow($sn, 3) * $c * pow($ok, 3));
        /** @var double $t16 */
        $t16 = 1.0 * (
                5.0
                + 6.0 * $eta
                + 28.0 * pow($t, 2)
                - 3.0 * pow($eta, 2)
                + 8.0 * pow($t, 2) * $eta
                + 24.0 * pow($t, 4)
                - 4.0 * pow($eta, 3)
                + 4.0 * pow($t, 2) * pow($eta, 2)
                + 24.0 * pow($t, 2) * pow($eta, 3)
            ) / (120.0 * pow($sn, 5) * $c * pow($ok, 5));
        /** @var double $t17 */
        $t17 = 1.0 * (61.0 + 662.0 * pow($t, 2) + 1320.0 * pow($t, 4) + 720.0 * pow($t, 6)) / (5040.0 * pow($sn, 7) * $c * pow($ok, 7));
        /** @var double $dlam */
        $dlam = $de * $t14 - pow($de, 3) * $t15 + pow($de, 5) * $t16 - pow($de, 7) * $t17;
        $lon = $olam + $dlam;
        $lon *= RAD2DEG;
        $lat *= RAD2DEG;
    }

    /**
     * Funkcja do konwersji wspórzędnych lat/lon WGS 84 na X/Y UTM
     *
     * @param int &$utmXZone
     * @param string &$utmYZone
     * @param double &$easting
     * @param double &$northing
     * @param double $lat
     * @param double $lon
     * @void
     */
    public function LatLonToUtmWGS84(&$utmXZone, &$utmYZone, &$easting, &$northing, $lat, $lon)
    {
        $this->LatLonToUtm(6378137.0, 1 / 298.257223563, $utmXZone, $utmYZone, $easting, $northing, $lat, $lon);
    }

    /**
     * Funkcja do konwersji wspórzędnych lat/lon WGS 84 na X/Y 1992 lub 2000
     *
     * @param double &$easting
     * @param double &$northing
     * @param double $lat
     * @param double $lon
     * @param int $proj
     * @void
     */
    public function LatLonToPUWGWGS84(&$easting, &$northing, $lat, $lon, $proj)
    {
        //$proj = 1 - dla odwzorowania kartograficznego 1992, każda inna wartość dla 2000
        $this->LatLonToPUWG(6378137.0, 1 / 298.257223563, $easting, $northing, $lat, $lon, $proj);
    }

    /**
     * Funkcja do konwersji wspórzędnych X/Y UTM na lat/lon elipsoidalne WGS 84
     *
     * @param int $utmXZone
     * @param string $utmYZone One character
     * @param double $easting
     * @param double $northing
     * @param double &$lat
     * @param double &$lon
     * @void
     */
    public function UtmToLatLonWGS84($utmXZone, $utmYZone, $easting, $northing, &$lat, &$lon)
    {
        $this->UtmToLatLon(6378137.0, 1 / 298.257223563, $utmXZone, $utmYZone, $easting, $northing, $lat, $lon);
    }

    /**
     * Funkcja do konwersji wspórzędnych X/Y 1992 lub 2000 na lat/lon elipsoidalne WGS 84
     *
     * @param double $easting
     * @param double $northing
     * @param int $proj
     * @param double &$lat
     * @param double &$lon
     * @void
     */
    public function PUWGToLatLonWGS84($easting, $northing, $proj, &$lat, &$lon)
    {
        $this->PUWGToLatLon(6378137.0, 1 / 298.257223563, $easting, $northing, $proj, $lat, $lon);
    }

}
