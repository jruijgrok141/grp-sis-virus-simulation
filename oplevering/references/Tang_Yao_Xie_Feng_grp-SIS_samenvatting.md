# Samenvatting: Tang, Yao, Xie, Feng (grp-SIS)

## Paper
**Titel:** *SIS Epidemic Modelling on Homogeneous Networked System: General Recovering Process and Mean-Field Perspective*  
**Auteurs:** Jiexi Tang, Yichao Yao, Meiling Xie, Minyu Feng (2025, arXiv:2505.12290)

## Kernidee
De auteurs breiden het klassieke SIS-model uit naar een **general recovering process SIS (grp-SIS)** model. In plaats van een geheugenloze (exponentiele) hersteltijd, laat dit model **willekeurige hersteltijdverdelingen** toe. Daardoor kan men realistischer niet-Markoviaanse dynamiek modelleren, zoals heavy-tailed of lognormale herstelpatronen.

## Wat is nieuw?
- Het model beschrijft herstel met een algemene dichtheid \(w(t)\) en overlevingsfunctie \(F_0(t)=P(W>t)\).
- Omdat herstel niet geheugenloos is, houden de auteurs expliciet de **infectieleeftijd** (hoe lang een node al besmet is) bij.
- Ze leiden hiervoor een mean-field vergelijking af op een homogeen netwerk met gemiddelde graad \(\langle k \rangle\), inclusief een integraalterm die de volledige herstelgeschiedenis samenvat.

## Belangrijkste theoretische resultaten
- Voor de stationaire toestand (niet-absorberend) leiden ze een algemene vorm af voor de verdeling van infectieleeftijden:
  \[
  f_{T(\infty)}(\tau)=\frac{F_0(\tau)}{\mathbb{E}[W]}.
  \]
- Ze tonen dat de gemiddelde infectieleeftijd in stationaire toestand afhangt van hogere momenten van de hersteltijd:
  \[
  \mathbb{E}[T(\infty)] = \frac{\mathbb{E}[W^2]}{2\mathbb{E}[W]}.
  \]
  Dit impliceert dat zware staarten (grote \(\mathbb{E}[W^2]\)) de gemiddelde infectieduur in de populatie sterk kunnen vergroten.
- De effectieve transmissiedrempel blijft in hun homogene mean-field afleiding:
  \[
  \tau_c=\frac{1}{\langle k \rangle}, \quad \tau=\beta\,\mathbb{E}[W].
  \]
  Dus de drempelvorm blijft gelijk aan klassiek SIS, ondanks algemene herstelverdelingen.

## Simulatie en interpretatie
- Het paper beschrijft een node-centric event-driven simulatieprocedure voor grp-SIS.
- De kernboodschap is dat de **vorm van de hersteltijdverdeling** vooral transiënten, infectieleeftijdsstatistieken en stationaire karakteristieken beïnvloedt, terwijl de basisdrempel in dit homogene mean-field kader structureel gelijk blijft.

## Praktische relevantie
Het model is bruikbaar wanneer herstelprocessen niet goed door een exponentiele verdeling worden benaderd (bijv. heterogene patiënttrajecten, bursty activiteitspatronen, of variabele “hersteltijd” in cyber-epidemieën). Daarmee biedt het een realistischer alternatief voor klassiek SIS zonder de analysemogelijkheden volledig te verliezen.

## Korte kritische noot
De resultaten zijn afgeleid voor een **homogene netwerk-aanname** en mean-field perspectief. Voor sterk heterogene of temporele netwerken kunnen drempel- en prevalentie-eigenschappen afwijken; de auteurs noemen dit ook als richting voor toekomstig werk.

---

## Kritische vergelijking met eigen resultaten (`report.pdf`)

Dit blok vergelijkt Tang et al. (grp-SIS, homogeen mean-field) met het eigen project *Testing a spectral epidemic threshold for SIS on networks under heavy-tailed recovery* (NetLogo ABM op één vaste ER-graaf, \(N\approx 2000\), \(\mathbb{E}[W]=5\), vergelijkbare herstelwetten).

### Waar theorie en simulatie elkaar aanvullen

- **Vorm van herstel vs. “drempel” in ruime zin.** Tang et al. benadrukken dat een algemene herstelverdeling vooral **niet-Markoviaanse structuur** (infectieleeftijd, hogere momenten van \(W\), stationaire infectieleeftijds-PDF) meeneemt. Jouw resultaten sluiten daar inhoudelijk bij aan: het **verschil tussen herstelwetten** is vooral zichtbaar in **extinctietijden en dynamiek** (median `bs-out-final-tick`; power-law langer dan exponentieel/lognormaal), terwijl de **verschuiving van \(\hat{\beta}_{\text{surv }50}\)** tussen wetten klein is vergeleken met de kloof tot de spectrale referentie.
- **Geen eenduidige “heavy tail maakt alles makkelijker persistent”.** Jouw hypothese (makkelijker persistentie bij zware staarten bij vaste \(\mathbb{E}[W]\)) wordt genuanceerd: lognormaal valt op het raster samen met exponentieel (gecensureerd bovenaan het \(\beta\)-grid); power-law kruist 50% survival iets **lager** dan de andere twee. Dat past bij het idee dat “heavy-tailed recovery” geen enkelvoudige hefboom is, maar per gekozen familie en implementatie anders werkt — vergelijkbaar met hoe Tang et al. de verdeling vooral in **statistische grootheden** (zoals \(\mathbb{E}[T(\infty)]\)) laten doorschemeren, niet in een grove drempelregel voor jouw ABM-outcome.

### Waar je resultaten het Tang-paper relativeren of scherp stellen

- **Andere referentie-drempel.** In Tang et al. verschijnt in het homogene mean-field kader \(\tau_c = 1/\langle k\rangle\) (structuur via gemiddelde graad). Jij gebruikt bewust de **spectrale** schaal \(\tau_{\text{pred}} = 1/\lambda_{\max}\) (en \(\beta_{\text{pred}} = 1/(\lambda_{\max}\mathbb{E}[W])\)). Dat is een andere — doorgaans scherpere — mean-field familie (NIMFA-achtig) dan louter \(1/\langle k\rangle\). Een letterlijke vergelijking “past Tang’s \(\tau_c\) op mijn \(\hat{\tau}_{50}\)?” vereist dus expliciet **welke theoretische lijn** je als tegenhanger kiest; jouw \(\tau_{\text{het}}=\langle k\rangle/\langle k^2\rangle\) verschuift slechts ~3% t.o.v. \(\tau_{\text{pred}}\) en verklaart de grote ABM-kloof niet.
- **Deterministisch endemic vs. stochastische extinctie op eindige horizon.** Tang et al. werken uit naar **stationaire** grootheden in een niet-absorberend scenario (eindige \(\rho_I^\infty\), verdeling van infectieleeftijd). Jouw **50%-survivaldrempel** meet of een **stochastisch** proces op **max-ticks** nog besmetting heeft — dat is conceptueel nabij “persistentie”, maar niet hetzelfde als het eindige-\(N\) absorberende SIS en het \(\mathbb{P}(\text{extinctie})\)-landschap. De grote overschotfactor \(\hat{\tau}_{50}/\tau_{\text{pred}} \approx 1{,}58\)–\(1{,}72\) weerspiegelt dus vooral **eindige \(N\)**, **discrete tijd**, correlaties buiten mean-field, en de gekozen survival-definitie — niet per se een tegenspraak met Tang’s **ODE/mean-field** voor \(\rho_I(t;\tau)\).
- **Discrete tijd vs. continue formulering.** Tang’s PDE-structuur is continue in de tijd; jouw model gebruikt **per-tick infectiekans** \(\beta\) en \(\tau = \beta\,\mathbb{E}[W]\) als discrete analogon. Kleine verschuivingen in effectieve drempels zijn daarmee verwacht, ook al is de grp-SIS idee hetzelfde.
- **Micro vs. mean-field traject (experiment 08).** Jouw Figuur 5 laat zien dat de **in-model homogeneous mean-field curve** en de microscopische grp-prevalentie dicht op elkaar liggen zolang extinctie niet ingrijpt — dat ondersteunt de interpretatie in de discussie: de **kloof in Tabel 3** zit in **waar stochastische survival waarschijnlijk wordt**, niet in een fout tussen “Tang-theorie” en de ODE die je parallel simuleert.

### Kort oordeel

Tang et al. leveren een **netjes afgebakende mean-field theorie** voor grp-SIS op homogene netwerken, waarin de **basisschaal** \(\tau_c = 1/\langle k\rangle\) en \(\tau = \beta\mathbb{E}[W]\) optreedt en de **herstelvorm** vooral **distributie van infectieleeftijd en transiënten** kleurt. Jouw werk toont in een **realistischer beslissings-setting** (één grote ER-graaf, discrete ABM, spectrale referentie, eindige horizon) dat **praktische drempels** en **herstelwet-gevoeligheid** vooral door **stochasticiteit, topologie-eigenwaarden en meetdefinitie** worden gedomineerd; daarmee **kwalificeer** je de vereenvoudigde hypothese dat zware staarten op zichzelf de spectrale lijn sterk “naar beneden trekken”, terwijl je tegelijk Tang’s bredere punt bevestigt dat **de vorm van \(W\) dynamiek en extinctie** sterk kan beïnvloeden zonder dat dat automatisch samenvalt met een scherpe, eenduidige verschuiving van een mean-field drempel.
