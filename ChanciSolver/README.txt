Los archivos se corren con una ruta relativa a la ubicaci�n del ejecutable. Es decir que si el ejecutable est� en la carpeta build, dentro de esta carpeta  (Que espero que s�) el archivo deber�a correrse como ./ChanciSolver.exe ../"Nombre del archivo.txt" Suponiendo que est� en la carpeta ChanciSolver.

Formato del archivo (REVISAR 1DWells.txt de ejemplo): TODO EST� EN SISTEMA INTERNACIONAL!

Malla
Dimensiones X Y Z
Grosores
X X X X (m)
Y Y (m)
Z (m)
TOPS
X X X X Y Y Y Y (m)

ROCK
COMPRESSIBILITY VALOR (1/Pa)
REFERENCE_PRESSURE VALOR (Pa)
POROSITY
VALORES POROSIDADES (-)
ABSOLUTE_PERMEABILITY
X1 Y1 Z1 X2 Y2 Z2 ... EN TOTAL 3*X*Y*Z (m2)

FLUIDS # (De fluidos)
TYPE Water/Gas/Oil
VOLUME_FACTOR #_MEDIDAS
PRESI�N 1 VOLUME_FACTOR 1 - presiones (Pa) vol_factor (m3rc/m3sc) o eso creo
PRESION 2 VOLUME_FACTOR 2
...(La cantidad de medidas)
VISCOSITY #_MEDIDAS
PRESI�N 1 VISCOSIDAD 1 - presiones (Pa) visc (Pa.s)
PRESION 2 VISCOSIDAD 2
...(La cantidad de medidas)
STANDARD_CONDITIONS_DENSITY VALOR
INITIAL_PRESSURE
X*Y*Z Presiones (Pa)
INITIAL_SATURATION
X*Y*Z Saturaciones (-)
PRINCIPAL 1/0 (1 o 0) Siempre debe haber uno en 1

TIME_DELTA VALOR (s)
SIMULATION_TIME VALOR (s)

WELLS 1
TYPE Producer/Injector
RADIUS VALOR (m)
NUMBER_OF_PERFORATIONS VALOR
POSITION X Y Z 
SKIN_FACTOR VALOR (-)

OPERATIVE_CONDITIONS
WELL 1
TYPE FLOW/PRESSURE/SHUT (Shut a�n no se implementa :v)
VALUE VALOR (m3/s)
NEXT_TIME VALOR (s) 

WELL 1
TYPE PRESSURE
VALUE VALOR
NEXT_TIME 999999999999 Para acabar la simulaci�n hay que poner un valor mayor que el tiempo de simulaci�n (No cambia), sigue con el constraint que llevaba hasta que acaba.
