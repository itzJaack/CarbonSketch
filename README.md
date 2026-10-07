# CarbonSketch

Laboratorio molecolare in italiano: ricerca per nome, formula molecolare e massa molare, struttura condensata, topologica, Lewis e modello 3D ruotabile. Layout adattabile a computer e telefoni; viste singole, ingrandimento, numerazione della catena principale, copia e download delle formule.

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
- Gli altri nomi vengono cercati tramite PubChem PUG REST. Sono disponibili alias italiani per diversi nomi comuni; per nomi non riconosciuti in italiano può essere necessario il nome inglese. Il nome sistematico restituito permette di controllare l’identità. Le risposte con più identificativi vengono rifiutate come ambigue.
- Le formule condensate lineari vengono generate soltanto dove la connettività è descrivibile correttamente. Per strutture con anelli complessi si mostra una formula semisviluppata con gruppi di atomi e collegamenti espliciti. Le forme aromatiche sono rappresentazioni di Kekulé, non legami alternati fisicamente localizzati.
- Lewis mostra tutti gli idrogeni e le coppie elettroniche libere per elementi/valenze supportati. Gli altri casi riportano una limitazione esplicita. Stereochimica e stereoisomeri richiedono nomi univoci e la verifica della struttura restituita.
- Se un conformero 3D non esiste, la dicitura distingue la geometria illustrativa dal conformero PubChem. WebGL è necessario per la visualizzazione interattiva.

## File

`dist/index.html` e `styles.css`: interfaccia; `app.js`: ricerca e interazioni; `chemistry.js`: nomenclatura, dati e formule; `renderers.js`: disegni SVG e 3D; `dist/vendor/3Dmol-min.js`: motore 3D incluso localmente. `serve.mjs` serve soltanto per l’anteprima locale.

Il sito statico è pubblicabile usando il contenuto di `dist/`. La configurazione Sites è in `.openai/hosting.json`. Le ricerche recenti sono conservate nel browser del visitatore. La ricerca di nuovi composti trasmette il nome cercato a PubChem; non sono richieste chiavi API o credenziali.

## Fonti

- PubChem: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest — dati di riferimento inclusi recuperati il 7 ottobre 2026.
- 3Dmol.js: https://3dmol.org/ e https://3dmol.org/doc/GLViewer.html — distribuzione ufficiale inclusa in `dist/vendor/`.
- Font DM Sans e Manrope da Google Fonts; in assenza di connessione vengono usati i caratteri disponibili sul dispositivo.
