# Kanoniczna struktura drugorzędowa RNA

Pliki `bpseq.py` oraz `pseudoknot.py` to główna część zadania. Ich uruchomienie
(np. `python3 bpseq.py`) spowoduje wypisanie listy testów, których kod obecnie
nie przechodzi np. `test_1l2x.json` dla `bpseq.py`.

# Krótkie omówienie kodu

## `bpseq.py`

- W kodzie mamy dwie główne klasy: `BPSEQ` i `DotBracket`
- Testy polegają na tym, że:
  - Program parsuje pliki JSON wygenerowane przez DSSR
  - Na podstawie wpisów `.dbn.all_chain.bseq` i `.dbn.all_chains.sstr` tworzona
    jest instancja `DotBracket`
  - Obiekt typu `DotBracket` konwertowany jest na obiekt typu `BPSEQ`
  - Program sprawdza czy jest to poprawny `BPSEQ`

## `pseudoknot.py`

- W kodzie mamy kilka głównych klas: `Strand`, `Hairpin`, `Stem`, `Pseudoknot`
  i `BPSEQ`
- Klasa `BPSEQ` została rozszerzona względem poprzedniej, bo posiada metody
  `stems()`, `hairpins()`, `pseudoknots()`, których celem jest odnalezienie
  odpowiednich motywów (reprezentowanych przez odpowiednie klasy)
- Testy polegają na tym, że:
  - Program wczytuje plik BPSEQ i parsuje go do postaci obiektu typu `BPSEQ`
  - Przy użyciu metod `stems()`, `hairpins()` i `pseudoknots()` generowane są
    listy motywów
  - Program sprawdza czy jest to zrobione poprawnie

# Zadania

## `bpseq.py`

- [2 pkt] Zaimplementuj metodę `DotBracket.from_string()`, tak aby poprawnie
  odczytała indeksy sparowanych nukleotydów
- [1 pkt] Zaimplementuj metodę `DotBracket.to_bpseq()`, tak aby poprawnie
  przygotowała listę `entries` (trójek, które wykorzystuje BPSEQ)

## `pseudoknot.py`

- [1 pkt] Zaimplementuj metodę `Stem.forms_pseudoknot_with()` do określania czy
  dwa stemy są w relacji typu pseudowęzeł
- [2 pkt] Zaimplementuj metodę `BPSEQ.stems()` do znajdowania stemów
- [2 pkt] Zaimplementuj metodę `BPSEQ.hairpins()` do znajdowania motywów _hairpin_

## Testy

Testy obecne w plikach stanowią punkt odniesienia. Poprawne implementacje ww.
metod dadzą wynik 0 błędów.
