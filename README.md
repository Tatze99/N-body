# N-body

## Beschreibung der Header.cpp 
- *set_satellite_old*: Erstellt ausgehend von der aktuellen Position der Erde eine Sonde, dabei kann die Geschwindigkeit in xy-Richtung und in z-Richtung variiert werden.
- *set_satellite*: Erstellt ausgehend von der aktuellen Position der Erde eine Sonde, wobei die Geschwindigkeit des Satelliten übergeben wird, die er zum Zeitpunkt t=0 bräuchte (berücksichtigt die Geschwindigkeitsänderung der Erde auf der elliptischen Bahn) 
- *crash_check*: Berechnet die Abstände der Sonde zu allen Planeten und überprüft, ob der Abstand kleiner als der Planetenradius ist.
- *set_values*: Wertet einen eingelesenen String schreibt diesen in entsprechende Variablen hinein. Als Trennzeichen werden Semicolons (ohne Leerzeichen) verwendet.
- *initialize_objects*: Löscht alte Objekte und erstellt neue Vektoren für Ort, Geschwindigkeit, Masse und Radius und weißt mithilfe von *set_values* Werte zu.
- *rk5-step*: Implementierung des Cash-Karp Verfahrens mithilfe zwei rk4-Integratoren, die in 4. und 5. Ordnung konvergieren. Ausgegeben werden die Differenzen für beide Integratoren zur Bestimmung der adaptiven Schrittweite. 
- *sat_driver*: Führt die Integration für eine bestimmte Zahl von Sonden durch, löscht Satelliten, die *crash_check* erfüllen. Berechnung von geeigneten Startwerten, für die Sonden auf der Bahn des Zielobjekts landen. 

## Beschreibung der Integrators.cpp 
- *all_from_target*
- *initialize_satellites*: Initialisiert Satelliten relativ zur Erdposition (10, 100, 1000 Stück)
- *integrator*: Integriert die Bewegung der Planten + Sonde für einen bestimmten Zeitraum mithilfe adaptiver Schrittweitensteuerung ohne die Werte in eine Datei zu schreiben. 
- *driver*: Integriert die Bewegung der Planten + Sonde, erstellt eine Datei, wo die Daten gespeichert werden 
- *check_for_boundaries*: Sucht nach minimaler und maximaler Geschwindigkeit, damit die Sonden eine Plantenbahn erreichen 
- *calc_sat*: Führt *check_for_boundaries* aus und berechnet für die das gefundene Geschwindigkeitsintervall die Trajektorien der Satelliten 
- *calc_angle*: Berechnet den Startzeitpunkt, damit eine Sonde tatsächlich einen Vorbeiflug an einem Planeten durchführen kann. Dafür wird zunächst eine Winkelabschätzung zur Bestimmung eines geeigneten Startzeitraums durchgeführt. Durch Variation von Geschwindigkeit wird der Abstand zum Planten minimiert. 
- *optimize_vsat*: Optimierung der Startgeschwindigkeit über ein 3-Punkte Problem unter der Annahme, dass der Zusammenhang zwischen Geschwindigkeit und minimalem Plantenabstand näherungsweise parabolisch ist. Mit drei verschiedenen Trajektorien wird das Minimum der Parabel gesucht, um eine optimale Startgeschwindigkeit zu bestimmen. 
- *optimize_initial_angle*: Verändert den Winkel in z-Richtung gegenüber der Startgeschwindigkeit der Erde, um sich dem Planeten für einen Swing by noch weiter anzunähern.

## Beschreibung der Programmteile 
- *fwd*: Integration der Planeten mit Hilfe des Vorwärts-Euler Verfahrens. Sehr ungenau und langsam.
- *rk4*: Integration der Planeten mit adaptiver Schrittweitensteuerung und einer rk4-Routine 
- *lf*: Leap-Frog Verfahren. Ungenau, da keine adaptive Schrittweitensteuerung. 
- *sat*: Ausführung von *calc_sat*
- *angle*: Ausführung von *Calc_angle*, Übergabe von Umlaufzeiten, Hohmann-Transferzeit und minimaler Geschwindigkeit der Sonde für einen bestimmten Planeten. 
- *calc_t*: Integration der PLnaten ausgehend von "Input_tend.csv" Datei, wo die Startdaten für den New-Horizons Start stehen.
- *Swing*: Führt *optmize_initial_angle* und *optimize_vsat* aus und integriert anschließend mit den optimierten Startparametern die Trajektorie.
- *Swing-by*: Berechne eine Trajektorie für einen Swing-by für verschiedene Zielplaneten. Die Geschwindigkeiten wurden durch die Routine *angle* ermittelt.
- *NewHorizon*: Berechne die Trajektorie der Raumsonde "New Horizon" (originale Startdaten, dann Optimierung)

### Beschreibung der Input.csv
- enthält folgende Daten: (x,y,z,vx,vy,vz,m)
- Input2.csv: (x,y,z,vx,vy,vz,m,r) enthält Planeten und Monde
- Reihenfolge der Inputobjekte: Sonne, Merkur, Venus, Erde, Mars, Jupiter, Saturn, Uranus, Neptun, Pluto, Sonde, Mond (relativ zur Erde)

### Beschreibung der Orbits.csv 
- Enthält folgende Daten: (a,e,b, a(1-e), a(1+e), siderische Umlaufzeit, Hohmann-Transferzeit, minimale Startgeschwindigkeit der Sonde)

