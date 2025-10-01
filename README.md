# GPS Satellite Position

**EN:**  
A simple Haskell program that computes the position of a GPS satellite from broadcast ephemeris data and a given GPS time.  
The algorithm is based on the official specification **IS-GPS-200H**.  
The program outputs satellite coordinates in the **ECEF (Earth-Centered, Earth-Fixed)** frame.

**PL:**  
Prosty program w Haskellu obliczający pozycję satelity GPS na podstawie efemeryd nadawanych przez satelitę oraz zadanego czasu GPS.  
Algorytm opracowano na podstawie specyfikacji **IS-GPS-200H**.  
Program wyznacza współrzędne satelity w układzie **ECEF (Earth-Centered, Earth-Fixed)**.

---

## Example output / Przykładowy wynik

```
Pozycja ECEF satelity [m]:
X =   15678943.123456789
Y =  -20456789.987654321
Z =    1345678.543210987
```

## Notes / Uwagi

- Input time must be given in **GPS time** (not UTC).  
- The type `UTCTime` is used only as a container for date arithmetic.  
- Leap seconds are **not** handled here – the user must provide GPS time directly.  
- All code is contained in a single file: `GpsSatelitePosition.hs`.  
```

---

## Data flow / Przepływ danych

GPS Ephemeris + GPS Time  
            │
            ▼
   Orbital computations
 (Kepler, harmonic corrections,  
 inclination, ascending node)    
            │
            ▼
 Transformation to ECEF frame 
            │
            ▼
 Satellite position (X, Y, Z)  
