# turbofan_prop
La cartella contiene scripts di MATLAB che analizzano un ciclo turbogas di un motore turbofan ed i suoi rendimenti, oltre che fornire la possibilità di realizzare grafici relativi

Per approfondire il tema abbiamo elaborato alcune function di MatLab che, 
dati una serie di input, permettono per ogni tipo di propulsore turbofan (a 
flussi separati o associati, con o senza post bruciatore) di calcolarne prestazioni 
e rendimenti. Per non eccedere di zelo ci siamo 
limitati al caso di volo subsonico, in quanto per concludere il computo di volo 
supersonico avremmo dovuto considerare tutti gli effetti dovuti alla 
comprimibilità. Questo avrebbe reso necessario inserire tra i parametri tutte le 
caratteristiche geometriche della presa dinamica, il che abbiamo 
concordato essere eccessivo e reo di sviare il focus dell’analisi.

Le function sono: 
	- TRBFN_SEP.m (per flussi separati, senza possibilità di post bruciatore)
	- TRBFN_AS.m (per flussi associati, in cui è possibile attivare il post 
		bruciatore a patto di inserirne temperatura massima e rendimento 
		pneumatico)
	- Piano_Ts.m (in cui, dati due vettori di tutte le pressioni e tutte le 
		temperature è possibile disegnare il ciclo corrispondente nel piano 
		temperatura/entropia).

Per TRBFN_AS.m abbiamo fatto in modo che la soluzione dipendesse da un 
𝐵𝑃𝑅 fissato (con la necessità di ricavare il 𝛽𝐹 corrispondente), poiché ci è 
sembrato un vincolo di progetto più realistico del rapporto di compressione del 
Fan.

Gli script che rendono possibile ricavare i grafici sono: 
	- plotrendimenti.m (per il motore senza postcombustore)
	- plotrendimentiAB.m (per il motore con postbruciatore).
