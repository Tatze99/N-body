# N-body

## Beschreibung der Header.cpp 
- *set_satellite_old*: Erstellt ausgehend von der aktuellen Position der Erde eine Sonde, dabei kann die Geschwindigkeit in xy-Richtung und in z-Richtung variiert werden.
- *set_satellite*: Erstellt ausgehend von der aktuellen Position der Erde eine Sonde, wobei die Geschwindigkeit des Satelliten übergeben wird, die er zum Zeitpunkt t=0 bräuchte (berücksichtigt die Geschwindigkeitsänderung der Erde auf der elliptischen Bahn) 
- *crash_check*:
Input.csv:(x,y,z,vx,vy,vz,m)
Input2.csv:(x,y,z,vx,vy,vz,m,r) enthält Planeten und Monde
Rheienfolge der Inputobjekte:
1.Sonne
2.Merkur
3.Venus
4.Erde
5.Mars
6.Jupiter
7.Saturn
8.Uranus
9.Neptun
10.Pluto
11.Mond			!Position bezieht sich auf den zugehörigen Planeten!
12.Ganymed(Jupiter)	!Position bezieht sich auf den zugehörigen Planeten!
13.Io(Jupiter)		!Position bezieht sich auf den zugehörigen Planeten!
14.Europa(Jupiter)	!Position bezieht sich auf den zugehörigen Planeten!
15.Titania(Uranus)	!Position bezieht sich auf den zugehörigen Planeten!
