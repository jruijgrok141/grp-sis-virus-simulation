# Overzicht voorgestelde verbeteringen — Review Jan Ruijgrok (OU)

Bron: `Review_Jan_Ruijgrok_OU.pdf` — peer review op het voorstel *Testing a spectral epidemic threshold for SIS on networks under heavy-tailed recovery*.

## Algemene punten om te verbeteren (Part 1)

1. **Precisie van begrippen en parameters**  
   De review merkt vaag taalgebruik aan, onder meer rond “typically stays present” en “heavy-tailed”. Voor wetenschappelijke strengheid moeten exacte parameters en definities in het voorstel (of het eindrapport) staan, zodat meting en replicatie betrouwbaar zijn.

2. **Initieel besmettingsaandeel**  
   Het is onduidelijk welk percentage agenten aanvankelijk geïnfecteerd is en welk effect die keuze op de dynamiek kan hebben. Dat verdient een expliciete keuze en uitleg.

## Aanbevelingen voor het eindrapport (bij hergebruik van de proposal-tekst)

1. **Numerieke grenzen / operationalisatie**  
   Vervang kwalitatieve formuleringen zoals “typically stays present” door een **kwantitatieve definitie** (cut-offs, regels).

2. **Visuele controle van de recovery-verdelingen**  
   Voeg **PDF-plots** (probability density functions) toe voor de drie recovery-verdelingen om aan te tonen dat de **gemiddelden gelijk** zijn terwijl de **spreiding / variantie** verschilt.

3. **Risk analysis in de discussie**  
   De risk analysis kan in de **discussie** worden verweven, mits er wordt uitgelegd **waarom** de voorgestelde mitigaties wel of niet werken.

## Aanvullende verbeteringen (uit andere vragen in het formulier)

- **Operationalisatie**: de onderzoeksvragen sluiten aan op meetbare simulatie-uitkomsten, maar de **algemene operationalisatie** kan scherper (aansluitend bij de punten over vage termen en drempels).
- **Detailniveau / reproduceerbaarheid**: het voorstel is sterk op software en theorie, maar **minder specifiek op simulatie-instellingen en modelparameters**; het eindrapport zou **meer concrete settings** moeten bevatten voor goede reproduceerbaarheid.
- **Kleine redactionele opmerking**: in het beoordeelde document wordt een typo “Pyton” genoemd (pagina 4).

## Korte samenvatting

De review benadrukt vooral **scherp definiëren en meten** (drempel/persistentie, heavy-tailed parameters, startconditie), **transparant maken van de recovery-vergelijking** (plots), en **meer implementatie- en parameterdetail** in het rapport, plus een **doordachte koppeling van risico’s en mitigatie** in de discussie.
