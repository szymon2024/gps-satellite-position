GPS SATELLITE POSITION

EN:  
A simple Haskell program that computes the position of a GPS satellite from broadcast ephemeris data and a given GPS time.  
The algorithm is based on the official specification IS-GPS-200H.  
The program outputs satellite coordinates in the ECEF (Earth-Centered, Earth-Fixed) frame.

PL:  
Prosty program w Haskellu obliczający pozycję satelity GPS na podstawie efemeryd nadawanych przez satelitę oraz zadanego czasu GPS.  
Algorytm opracowano na podstawie specyfikacji IS-GPS-200H.  
Program wyznacza współrzędne satelity w układzie ECEF (Earth-Centered, Earth-Fixed).


Example output / Przykładowy wynik

ECEF satellite position [m]:
X = 22121977.179040890
Y = 13231276.153471118
Z =  7440056.739046174


NOTES / UWAGI

EN:
- Input time must be given in GPS time (not UTC).
- The type UTCTime is used only as a container for date arithmetic.  
- Leap seconds are NOT handled here - the user must provide GPS time directly.  
- All code is contained in a single file: GpsSatellitePosition.hs.

PL:
- Czas wejściowy musi być podany w czasie GPS (nie UTC).
- Typ UTCTime służy jedynie jako kontener dla arytmetyki daty.
- Sekundy przestępne NIE są tutaj obsługiwane – użytkownik musi podać czas GPS bezpośrednio.
- Cały kod znajduje się w jednym pliku: GpsSatellitePosition.hs.


DATA FLOW / PRZEPŁYW DANYCH

GPS Ephemeris + GPS Time  
            |
            V
   Orbital computations  
 (Kepler, harmonic corrections,  
 inclination, ascending node)  
            |
            V
 Transformation to ECEF frame  
            |
            V
 Satellite position (X, Y, Z)  
