# Zadanie 1 - deduplikacja pachołków i najkrótsza otoczka

Program w C++ wczytuje listę wykrytych pachołków, scala odczyty leżące bliżej siebie niż zadana minimalna odległość (1.1), a następnie wyznacza otoczkę o najkrótszym obwodzie (1.2) na zbiorze po deduplikacji.

## Wymagania

- kompilator z obsługą **C++23** (testowane na GCC 13.3),
- CMake ≥ 3.16 (opcjonalnie).

## Budowanie i uruchomienie

CMake:

```bash
cmake -S . -B build
cmake --build build
./build/src/zad1 data/zad1/example_input.txt
```

Bez CMake:

```bash
g++ -std=c++23 -Wall -Wextra -O2 src/zad1.cpp -o zad1
./zad1 data/zad1/example_input.txt
```

Program przyjmuje dokładnie jeden argument: ścieżkę do pliku wejściowego. Przy błędzie wypisuje komunikat na `stderr` i kończy się kodem `1`.

## Format danych

**Wejście** (plik `.txt`):

```
<n>            liczba pachołków
<m>            minimalna dopuszczalna odległość między pachołkami
<id> <x> <y>   n linii: unikalne ID i współrzędne, oddzielone spacjami
```

**Wyjście:**

```
CONE <nowe_id> <x> <y> [<id składowych>]      dla każdego nowego (scalonego) pachołka
ENCLOSURE <obwód> [<id pachołków otoczki>]
```

## Rozwiązanie

### Struktura programu

| Funkcja | Rola |
|---|---|
| `readInput` | wczytanie i walidacja pliku (`getline` + `istringstream`, linia po linii) |
| `mergeAll` | deduplikacja (1.1) |
| `findClosestPair` | najbliższa para klastrów w bieżącym stanie |
| `assignIds` | nadanie ID po deduplikacji |
| `toCones` | konwersja klastrów na pachołki |
| `getBorderPath` | otoczka wypukła (1.2) |
| `calculate_perimeter` | obwód otoczki |

### 1.1 Deduplikacja

Zastosowany algorytm to **aglomeracyjna klasteryzacja hierarchiczna z centroid linkage** i progiem zatrzymania `m`:

1. Każdy pachołek jest na starcie osobnym klastrem.
2. Wyszukiwana jest para klastrów o najmniejszej odległości między środkami.
3. Jeśli ta odległość jest `≥ m`, algorytm się kończy. W przeciwnym razie para zostaje scalona i następuje powrót do kroku 2.

**Dlaczego najpierw najbliższa para.** Wynik zachłannego scalania „dowolnej pary bliżej niż `m`” zależy od kolejności przetwarzania. Dla punktów na prostej A=0, B=0.9, C=1.8 i m=1 scalenie A+B daje {0.45, 1.8}, a scalenie B+C daje {0, 1.35}. Wybór zawsze najbliższej pary usuwa tę zależność: najpierw scalane są odczyty, które najprawdopodobniej są duplikatem.

**Gwarancja na wyjściu.** Po zakończeniu żadne dwa pachołki nie leżą bliżej niż `m`. Każde scalenie zmniejsza liczbę klastrów o 1, więc pętla wykonuje się najwyżej `n − 1` razy.

**Uśrednianie ważone.** Klaster przechowuje sumy współrzędnych (`sumX`, `sumY`) i listę oryginalnych ID, a nie sam środek. Scalenie polega na dodaniu sum, a środek to `suma / liczba członków`, czyli średnia ze wszystkich oryginalnych odczytów. Uśrednianie samych środków byłoby błędne dla klastrów różnej wielkości: dla klastra {4,10} + 9 dałoby x = 5.175 zamiast poprawnego 5.067.

**Nadawanie ID.** Nowe ID nadawane są dopiero po zakończeniu deduplikacji, żeby scalenia pośrednie (np. 4+10 przed dołączeniem 9) nie zużywały numerów. Numeracja zaczyna się od `max(id) + 1`, a nie od `n + 1`, bo treść gwarantuje tylko unikalność ID, a nie ich ciągłość od 1. Scalone klastry są numerowane w kolejności najmniejszego ID składowej, więc wynik nie zależy od kolejności linii w pliku.

**Szczegóły implementacji:**
- porównywane są kwadraty odległości (`d² < m²`), co pozwala uniknąć `sqrt` w pętli i nie zmienia porządku,
- usuwanie scalonego klastra metodą *swap-and-pop* w O(1), bo kolejność klastrów w wektorze nie ma znaczenia,
- warunek scalenia to `d < m`: skoro `m` jest *minimalną dopuszczalną* odległością, to `d = m` jest dozwolone.

**Złożoność:** O(n³) (najwyżej n−1 scaleń, każde z wyszukiwaniem pary w O(n²)). Dla liczby pachołków spotykanej na torze to w zupełności wystarczy.

**Rozważane alternatywy:**
- *Single-linkage* (spójne składowe grafu „d < m”) nie zależy od kolejności, ale łańcuch odczytów może rozciągnąć się dowolnie daleko, a środki sąsiednich składowych nie muszą spełniać warunku `m`.
- *Kolejka priorytetowa* z leniwym usuwaniem nieaktualnych par obniża złożoność do O(n² log n). Wystarczy wrzucać do niej tylko pary o `d < m`.
- *Siatka przestrzenna* o boku `m`: duplikaty mogą leżeć tylko w sąsiednich komórkach 3×3, co ogranicza liczbę porównań.
- Algorytm *nearest-neighbor chain* nie nadaje się tutaj, bo działa tylko dla redukowalnych linkage'ów, a centroid linkage do nich nie należy.

### 1.2 Najkrótsza otoczka

**Najkrótsza otoczka to otoczka wypukła.** Każdy wielokąt ogradzający wszystkie punkty musi zawierać ich otoczkę wypukłą, a zbiór wypukły ma najmniejszy obwód spośród wszystkich krzywych zamkniętych, które go zawierają.

Do jej wyznaczenia użyto algorytmu **Andrew's monotone chain**, O(n log n):

1. Punkty są sortowane leksykograficznie (po `x`, przy remisie po `y`).
2. Budowana jest dolna otoczka od lewej do prawej, a potem górna od prawej do lewej. Za każdym razem odrzucane są punkty, w których ścieżka nie skręca w lewo.
3. Kierunek skrętu określa iloczyn wektorowy: `> 0` oznacza skręt w lewo, `< 0` w prawo, a `= 0` punkty współliniowe.

Obwód to suma długości boków, łącznie z bokiem domykającym (ostatni wierzchołek → pierwszy).

Otoczka jest liczona na zbiorze **po deduplikacji**: scalone pachołki zastępują swoje składowe.

## Edge case'y

| Przypadek | Obsługa |
|---|---|
| plik nie istnieje | błąd: `Unable to open a file ...` |
| `n` ujemne lub niecałkowite | błąd linii 1; `n` jest wczytywane do typu ze znakiem, bo `"-5"` wczytane do `size_t` zawinęłoby się do ogromnej liczby zamiast zgłosić błąd |
| `m ≤ 0` | błąd linii 2 (odległość musi być dodatnia) |
| za mało linii, brakujące pole, tekst zamiast liczby, nadmiarowe pole | błąd z numerem linii |
| spacje na końcu linii, `\r` z Windowsa | akceptowane (`std::ws`) |
| pusta linia | błąd formatu: wejście musi ściśle trzymać się specyfikacji |
| `d = m` dokładnie | pachołki nie są scalane |
| łańcuch odczytów A-B-C | rozstrzygany regułą najbliższej pary (patrz wyżej) |
| kaskadowe scalenia | środek liczony z sum, ID nadawane na końcu |
| ID nieciągłe (np. 5, 17, 100) | nowe ID od `max(id) + 1`, brak kolizji |
| remis odległości | wygrywa pierwsza znaleziona para (przy liczbach zmiennoprzecinkowych rzadkie) |
| ≤ 3 pachołki | wszystkie należą do otoczki |
| punkty współliniowe na krawędzi otoczki | pomijane na liście ID (warunek `≤ 0`); obwód bez zmian |
| wszystkie punkty na jednej prostej | otoczka zdegenerowana do odcinka, obwód „tam i z powrotem” |
| zduplikowane punkty | nie występują po deduplikacji, więc nie zakłócają monotone chain |

## Dowód poprawności

**Przykład z treści zadania** (`data/zad1/example_input.txt`, `m = 1.0`):

```
CONE 13 10.300 0.400 [2,8]
CONE 14 5.067 9.700 [4,9,10]
ENCLOSURE 39.162 [12,1,11,13,3,14,5]
```

Wynik zgadza się z oczekiwanym: te same nowe pachołki, ten sam obwód i ten sam zbiór siedmiu ID otoczki (kolejność na liście według treści nie jest istotna).

Przebieg deduplikacji krok po kroku:

| Krok | Najbliższa para | Odległość | Wynik |
|---|---|---|---|
| 1 | 4 - 10 | 0.583 | {4,10}, środek (4.85, 9.75) |
| 2 | 2 - 8 | 0.632 | {2,8}, środek (10.3, 0.4) |
| 3 | {4,10} - 9 | 0.667 | {4,9,10}, środek (5.067, 9.700) |
| 4 | - | > 1.0 | koniec |

Krok 3 pokazuje, dlaczego odległości trzeba liczyć od nowa po każdym scaleniu: odległości 0.667 nie było wśród oryginalnych par.

**Niezależność od kolejności danych:** ten sam plik z odwróconą kolejnością linii pachołków daje identyczne wyjście.

**Walidacja wejścia:** plik z linią `2 5` (brak współrzędnej `y`) kończy się komunikatem `Error: Line 4: Wrong number` i kodem wyjścia `1`.