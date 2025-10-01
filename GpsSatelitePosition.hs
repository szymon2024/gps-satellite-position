{-# LANGUAGE RecordWildCards #-}

{- EN: Determining the GPS satellite position from the GPS ephemeris and a given GPS time.
   EN: Based on IS-GPS-200H.
   PL: Wyznaczenie pozycji satelity GPS z efemeryd GPS i danego czasu GPS.
   PL: Opracowano na podstawie IS-GPS-200H.
-}

import Data.Time.Calendar  (Day(..), fromGregorian, fromGregorianValid)
import Data.Time.Clock     (UTCTime(..), DiffTime, NominalDiffTime, diffUTCTime)
import Data.Time.LocalTime (makeTimeOfDayValid, timeOfDayToTime)
import Data.Fixed          (Pico)
import Text.Printf         (printf)

-- | EN: GPS ephemeris (a subset of fields necessary to determine a satellite's position)
-- | PL: Efemerydy GPS (podzbiór pól niezbędnych do wyznaczenia pozycji satelity)
data Ephemeris = Ephemeris
  { crs      :: Double            -- ^ EN: orbital radius correction [m]
                                  -- ^ PL: poprawka harmoniczna promienia [m]
  , deltaN   :: Double            -- ^ EN: mean motion difference [rad/s]
                                  -- ^ PL: średnia różnica ruchu [rad/s]
  , m0       :: Double            -- ^ EN: mean anomaly at reference epoch [rad]
                                  -- ^ PL: anomalia średnia w chwili odniesienia [rad]
  , cuc      :: Double            -- ^ EN: latitude argument correction [rad]
                                  -- ^ PL: poprawka harmoniczna argumentu szerokości [rad]
  , e        :: Double            -- ^ EN: eccentricity
                                  -- ^ PL: ekscentryczność
  , cus      :: Double            -- ^ EN: latitude argument correction
                                  -- ^ PL: poprawka harmoniczna argumentu szerokości [rad]
  , sqrtA    :: Double            -- ^ EN: sqare root of semi-major axis [m^0.5]
                                  -- ^ PL: pierwiastek kwadratowy z półosi wielkiej [m^0.5]
  , toe      :: Double            -- ^ EN: time of ephemeris in GPS week [s]
                                  -- ^ PL: czas efemeryd w tygodniu GPS [s]
  , cic      :: Double            -- ^ EN: inclination correction [rad]
                                  -- ^ PL: poprawki harmoniczne inklinacji [rad]
  , omega0   :: Double            -- ^ EN: longitude of ascending node at the beginning of the week [rad]
                                  -- ^ PL: długość geograficzna węzła wstępującego na początku tygodnia [rad]
  , cis      :: Double            -- ^ EN: inclination correction [rad]
                                  -- ^ PL: poprawka harmoniczna inklinacji [rad]
  , i0       :: Double            -- ^ EN: inclination at reference epoch [rad]
                                  -- ^ PL: nachylenie orbity w epoce odniesienia [rad]
  , crc      :: Double            -- ^ EN: orbital radius corrcetion [m]
                                  -- ^ PL: poprawka harmoniczna promienia [m]
  , omega    :: Double            -- ^ EN: argument of perigee
                                  -- ^ PL: argument perygeum
  , omegaDot :: Double            -- ^ EN: rate of node's right ascension
                                  -- ^ PL: szybkość rektascensji węzła [rad/s]
  , iDot     :: Double            -- ^ EN: rate of inclination angle                                                    
  } deriving (Show)

data GPSTime = GPSTime { gpsDay     :: Day
                       , gpsDayTime :: DiffTime
                       }
  deriving (Eq, Ord, Show)
               
-- | EN: Constants
-- | PL: Stałe
mu, omegaE :: Double
mu     = 3.986005e14            -- EN: value of Earth's universal gravitational parameters in WGS84 [m^3/s^2]
                                -- PL: stała grawitacyjna Ziemi w WGS84 [m^3/s^2]
omegaE = 7.2921151467e-5        -- EN: Earth's rotational speed [rad/s]
                                -- PL: szybkość obrotowa Ziemi w WGS84 [rad/s]

-- | EN: Determining the GPS satellite position from the GPS ephemeris and for a given GPS time
-- | PL: Wyznaczenie pozycji satelity z efemeryd i dla danego czasu
gpsSatellitePosition         
    :: Double                             -- ^ EN: time in seconds of week [s]
    -> Ephemeris                          -- ^ EN: ephemeris
    -> (Double, Double, Double)           -- ^ EN: (X, Y, Z) in ECEF
gpsSatellitePosition  t Ephemeris{..} =
  let
    a  = (sqrtA)*(sqrtA)                  -- EN: semi-major axis
                                          -- PL: półoś wielka 
    n0 = sqrt(mu/(a*a*a))                 -- EN: computed mean motion [rad/sec]
                                          -- PL: średni ruch [rad/sec]
    n  = n0 + deltaN                      -- EN: corrected mean motion [rad/s]
                                          -- PL: skorygowany ruch średni [rad/s]
    tk = wrapWeekCrossover (t - toe)      -- EN: time from ephemeris reference epoch toe [s]
                                          -- PL: czas od epoki odniesienia efemeryd toe [s]
    mk = m0 + n*tk                        -- EN: mean anomaly for tk
                                          -- PL: średnia anomalia
    ek = keplerSolve mk e                 -- EN: eccentric anomaly [rad]
                                          -- PL: anomalia mimośrodowa [rad]
                                          
    vk = atan2 (sqrt (1 - e*e) * sin ek) (cos ek - e)  -- EN: true anomaly
                                                       -- PL: anomalia prawdziwa
    phik = vk + omega                                  -- EN: argument of latitude
                                                       -- PL: argument szerokości geograficznej
                                                        
             
    duk  = cus * sin (2*phik) + cuc * cos (2*phik)     -- EN: argument of latitude correction
                                                       -- PL: argument korekcji szerokości geograficznej
    drk  = crs * sin (2*phik) + crc * cos (2*phik)     -- EN: radius correction
                                                       -- PL: korekta promienia
    dik  = cis * sin (2*phik) + cic * cos (2*phik)     -- EN: inclination correction
                                                       -- PL: korekta inklinacji
             
    uk     = phik + duk                                -- EN: corrected argument of latitude
                                                       -- PL: skorygowany argument szerokości geograficznej
                                         
    rk     = a * (1 - e*cos ek) + drk                  -- EN: corrected radius
                                                       -- PL: skorygowany promień
    ik     = i0 + dik + iDot*tk                        -- EN: corrected inclination
                                                       -- PL: skorygowana inklinacja
             
    xk  = rk * cos uk                                  -- EN: positions in the orbital plane
    yk  = rk * sin uk                                  -- PL: pozycje w płaszczyźnie orbitalnej
             
    omegak = omega0 + (omegaDot - omegaE)*tk - omegaE*toe -- EN: corrected longitude of ascending node
                                                          -- PL: skorygowana długość geograficzna węzła wstępującego
    xk' = xk * cos omegak - yk * cos ik * sin omegak      -- EN: transformation to ECEF
    yk' = xk * sin omegak + yk * cos ik * cos omegak      -- PL: transformacja do ECEF
    zk' =                   yk * sin ik
  in (xk',yk',zk')


-- | EN: Iterative solution of Kepler's equation ek = m + e sin ek (Metoda Newtona-Raphsona)
-- | PL: Rozwiązanie iteracyjne równania Keplera ek = m + e sin ek (Metoda Newtona-Raphsona)
keplerSolve    
    :: Double  -- ^ EN: mean anomaly
    -> Double  -- ^ EN: eccentricity
    -> Double  -- ^ EN: eccentric anomaly [rad]
keplerSolve m e = go e0 0               
  where
    e0 = m + e * sin m
    go :: Double -> Int -> Double
    go eN k
      | k > 20 = eN
      | abs (eN' - eN) < 1e-12 = eN'
      | otherwise = go eN' (k+1)
          where    
            f     = eN - e * sin eN - m  
            fp    = 1 - e * cos eN       -- EN: fp derivative of the function f
                                         -- PL: fp pochodna funkcji f
            eN'   = eN - f/fp            -- EN: iterative formula
                                         -- PL: formuła iteracyjna

fullWeek, halfWeek :: Double
fullWeek = 604800.0  -- EN: number of seconds of full week
                     -- PL: liczba sekund całego tygodnia
halfWeek = 302400.0  -- EN: number of seconds of half week
                     -- PL: liczba sekund połowy tygodnia

-- | EN: Normalization to +/- half a week
-- | PL: Normalizacja do +/- pół tygodnia
wrapWeekCrossover
    :: Double      -- ^ EN: time difference of GPS week [s]
    -> Double      -- ^ EN: time difference of GPS half week [s]
wrapWeekCrossover dt
  | dt >  halfWeek = dt - fullWeek
  | dt < -halfWeek = dt + fullWeek
  | otherwise      = dt

-- | EN: Creates GPSTime from given numbers (interpreted as GPS time, not UTC!)
-- | PL: Tworzy GPSTime z podanych liczb (interpretowanych jako czas GPS, nie UTC!)
mkGPSTime :: Integer -> Int -> Int -> Int -> Int -> Pico -> Maybe GPSTime
mkGPSTime year month day h m s = do
  tod <- makeTimeOfDayValid h m s
  d   <- fromGregorianValid year month day
  return $ GPSTime d (timeOfDayToTime tod)

-- | EN: GPS epoch: 1980-01-06 00:00:00 UTC
-- | PL: Epoka GPS: 1980-01-06 00:00:00 UTC
gpsEpoch :: UTCTime
gpsEpoch = UTCTime (fromGregorian 1980 1 6) 0

-- | EN: Conversion GPSTime -> (week, sow).
-- |     Uses UTCTime only as a container for date arithmetic.
-- |     Input must already be GPS time (no leap second correction applied here).
-- | PL: Konwersja GPSTime -> (week, sow).
-- |     UTCTime używany jest wyłącznie jako kontener do obliczeń na datach.
-- |     Wejście musi być już w czasie GPS (bez korekcji sekund przestępnych).
gpsTimeToWeekSow :: GPSTime -> (Integer, NominalDiffTime)
gpsTimeToWeekSow (GPSTime d t) =
    let gpsTime' = UTCTime d t
        diffSec  = diffUTCTime gpsTime' gpsEpoch
        week = floor (diffSec / 604800)
        sow  = diffSec - fromIntegral week * 604800
    in (week, sow)
     
-- | EN: Ephemeris example
-- | PL: Przykład efemeryd
ephsExample :: Ephemeris
ephsExample = Ephemeris
          { crs = 134.75
          , deltaN = 4.12374319903e-9
          , m0 = -2.23037324521
          , cuc = 7.25500285625e-6
          , e = 1.61300709005e-2
          , cus = 8.00378620625e-6
          , sqrtA = 5153.71900558
          , toe = 424800.0
          , cic = 7.07805156708e-8
          , omega0 = -2.84371723002
          , cis = 2.6635825634e-7
          , i0 = 0.967821256634
          , crc = 226.75
          , omega = -1.2251476939
          , omegaDot = -7.73246494535e-9
          , iDot = 2.60367988229e-10
          }

{-
-- | EN: Expected result for ephemeris example
-- | PL: Oczekiwany wynik dla przykładu efemeryd
expectedResult :: (Double, Double, Double)
expectedResult = ( 22121977.179040890
                 , 13231276.153471118
                 , 7440056.739046174
                 )
-}

-- | EN: Calculate GPS satelite position for GPS ephemeris example and GPS time
-- | PL: Oblicza pozycję satelity GPS dla przykładu efemeryd GPS i czasu GPS
main :: IO ()
main = do
  let Just gpsT  = mkGPSTime 2024 03 07 21 59 30.0
      ephs       = ephsExample
      (_,sow)    = gpsTimeToWeekSow gpsT
      (xS,yS,zS) = gpsSatellitePosition (realToFrac sow) ephs
                
  putStrLn $ "Pozycja ECEF satelity [m]:"
  putStrLn $ "X = " ++ printf "%18.9f" xS
  putStrLn $ "Y = " ++ printf "%18.9f" yS
  putStrLn $ "Z = " ++ printf "%18.9f" zS

