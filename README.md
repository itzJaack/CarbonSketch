# CarbonSketch

Laboratorio molecolare in italiano: ricerca per nome, formula molecolare e massa molare, struttura condensata, condensata 2D con catena orizzontale e rami verticali, topologica, Lewis e modello 3D ruotabile. Layout adattabile a computer e telefoni; viste singole, ingrandimento, numerazione della catena principale, copia e download delle formule.

## Inserimento e disegno

Scrivi un nome in italiano o in inglese. L'interprete tollera spazi, trattini fra i locanti e trattini omessi fra numeri e sostituenti: `2-3-4-trimetilesano` viene letto come `2,3,4-trimetilesano`, senza ricorrere a una ricerca in inglese. La pagina mostra l'interpretazione adottata. Le correzioni riguardano la scrittura: locanti impossibili, valenze e quantità incoerenti dei sostituenti producono un errore.

`4,5dietil2,3,4,5tetraetilnonano` conserva tutti e sei i gruppi etilici della descrizione, anche quando più gruppi si trovano sullo stesso carbonio. Quando la catena indicata non è la più lunga, la pagina segnala che il nome non segue la scelta della catena principale IUPAC. Il nome ricostruito localmente non è una certificazione della nomenclatura completa.

Seleziona **Disegna condensata 2D** per posizionare gruppi e collegarli con legami singoli, doppi o tripli. Puoi partire dalla molecola attuale, modificare il disegno o leggere una formula come `CH3-CH2-CH(CH2-CH3)-CH2-CH(C(CH3)3)-CH2-CH2-CH3`. **Analizza disegno** controlla connessioni e valenze e genera tutte le rappresentazioni. I numeri di idrogeni scritti nei gruppi vengono rispettati e verificati. Le strutture non nominabili dal sottoinsieme locale sono chiamate "Molecola disegnata".

Il disegno analizzato viene conservato soltanto nel browser. La geometria 3D del disegno è una stima didattica, senza ottimizzazione energetica; il disegno non viene inviato a un servizio esterno per assegnare un nome.

## Pubblicazione

- Repository: https://github.com/itzJaack/CarbonSketch
- GitHub Pages: https://itzjaack.github.io/CarbonSketch/
- Sites: https://carbonsketch.itzjaack.chatgpt.site/

Il workflow `.github/workflows/pages.yml` esegue i controlli e pubblica il contenuto di `dist/` su GitHub Pages. La pubblicazione non include credenziali né una fase di compilazione. I percorsi degli asset sono relativi e funzionano anche sotto `/CarbonSketch/`.

## Avvio locale

Con Node.js installato, dalla cartella del progetto:

```
npm start
```

Apri http://127.0.0.1:4173. Non servono installazioni aggiuntive né una fase di compilazione. Per i controlli automatici: `npm test`.

## Nomi e correttezza chimica

- L’esempio `5-isobutil decano` usa il composto PubChem CID 21201047: C14H30, massa molare 198,39 g/mol, catena principale di dieci carboni e ramificazione CH2–CH(CH3)2 sul carbonio 5. Il sito include coordinate 2D e un conformero 3D PubChem.
- Anche etanolo, metanolo, benzene, acetone, acido acetico, acqua e caffeina includono dati PubChem nel progetto; questi esempi non richiedono richieste a un servizio esterno dopo il caricamento della pagina.
- Il riconoscimento locale comprende catene C1–C20, alcani ramificati con sostituenti comuni metilici, etilici, propilici, butilici, iso-, sec- e tert-, cicloalcani saturi e alcuni alcheni/alchini aciclici con un legame multiplo. I conformeri generati localmente sono sempre indicati come geometrie illustrative stimate, senza ottimizzazione energetica.
- I nomi fuori dal sottoinsieme locale vengono cercati tramite PubChem PUG REST. Sono disponibili alias italiani per diversi nomi comuni; i nomi complessi fuori dal sottoinsieme possono richiedere la forma inglese. Un risultato non trovato non prova l'inesistenza del composto. Il nome restituito permette di controllare l’identità. Le risposte con più identificativi vengono rifiutate come ambigue.
- Le formule condensate lineari vengono generate soltanto dove la connettività è descrivibile correttamente. Per strutture con anelli complessi si mostra una formula semisviluppata con gruppi di atomi e collegamenti espliciti. Le forme aromatiche sono rappresentazioni di Kekulé, non legami alternati fisicamente localizzati.
- Lewis mostra tutti gli idrogeni e le coppie elettroniche libere per elementi/valenze supportati. Gli altri casi riportano una limitazione esplicita. Stereochimica e stereoisomeri richiedono nomi univoci e la verifica della struttura restituita.
- Se un conformero 3D non esiste, la dicitura distingue la geometria illustrativa dal conformero PubChem. WebGL è necessario per la visualizzazione interattiva.

## File

`dist/index.html` e `styles.css`: interfaccia; `app.js`: ricerca e interazioni; `chemistry.js`: nomenclatura, dati e validazione; `condensed2d.js`: catena e rami in SVG; `condensed-editor.js`: disegno e lettura delle formule condensate; `renderers.js`: disegni SVG e 3D; `dist/vendor/3Dmol-min.js`: motore 3D incluso localmente. `serve.mjs` serve soltanto per l’anteprima locale.

Il sito statico è pubblicabile usando il contenuto di `dist/`. La configurazione Sites è in `.openai/hosting.json`. Le ricerche recenti sono conservate nel browser del visitatore. La ricerca di nuovi composti trasmette il nome cercato a PubChem; non sono richieste chiavi API o credenziali.

## Fonti

- PubChem: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest — dati di riferimento inclusi recuperati il 7 ottobre 2026.
- 3Dmol.js: https://3dmol.org/ e https://3dmol.org/doc/GLViewer.html — distribuzione ufficiale inclusa in `dist/vendor/`.
- Font DM Sans e Manrope da Google Fonts; in assenza di connessione vengono usati i caratteri disponibili sul dispositivo.
