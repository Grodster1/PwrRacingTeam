
/**
 * Zadanie 1 – deduplikacja pachołków i najkrótsza otoczka.
 *
 * Przebieg programu:
 *   1. readInput - wczytanie i walidacja pliku wejściowego,
 *   2. mergeAll - deduplikacja: klasteryzacja hierarchiczna
 *                 (centroid linkage) z progiem zatrzymania minDist,
 *   3. assignIds - nadanie ID: pojedyncze pachołki zachowują swoje ID,
 *                  scalone dostają nowe, począwszy od max(id) + 1,
 *   4. getBorderPath - otoczka wypukła (Andrew's monotone chain, O(n log n)),
 *   5. calculate_perimeter - obwód otoczki.
 *
 * Najkrótsza otoczka = otoczka wypukła: każdy wielokąt ogradzający wszystkie
 * punkty musi zawierać ich otoczkę wypukłą, a zbiór wypukły ma najmniejszy
 * obwód spośród krzywych zamkniętych, które go zawierają.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <expected>
#include <stdexcept>
#include <string>
#include <vector>
#include <optional>
#include <limits>
#include <algorithm>
#include <cmath>

// Pojedynczy pachołek: unikalne ID i pozycja na mapie.
struct Cone{
    int id;
    double x, y;

    friend std::ostream& operator<<(std::ostream& os, Cone const& c){
        return os << "Cone(id=" << c.id << ", pos=(" << c.x << ", " << c.y << "))";
    }

    // Porządek leksykograficzny (najpierw x, przy takim samym y),
    // wymagany do posortowania punktów.

    bool operator<(const Cone& c) const{
        return x < c.x || (x == c.x && y < c.y);
    }
};


// Dane z pliku wejściowego: deklarowana liczba punktów, próg m i lista pachołków.
struct InputData{
    std::size_t n;
    double minDist;
    std::vector<Cone> cones;

    friend std::ostream& operator <<(std::ostream& os, InputData const& data)
    {
        os << data.n << '\n'
           << data.minDist << '\n';
        for (const auto& cone : data.cones){
            os << cone << '\n';
        }
        return os;       
    }
};


/**
 * Grupa odczytów uznanych za jeden fizyczny pachołek.
 *
 * Zamiast środka przechowujemy sumy współrzędnych. Dzięki temu scalenie to
 * zwykłe dodanie sum, a środek (sum / liczba członków) jest automatycznie
 * średnią ważoną wszystkich oryginalnych punktów. Uśrednianie samych środków
 * dałoby zły wynik, gdy klastry mają różną liczbę członków.
 */

struct Cluster{
    double sumX, sumY;          
    std::vector<int> membersId;
    int id;

    double center_x() const {return sumX/membersId.size();};
    double center_y() const {return sumY/membersId.size();};

    friend std::ostream& operator <<(std::ostream& os, Cluster const& cluster)
    {
        os << "ID:" << cluster.id << " ";
        os << "Center: (" << cluster.center_x() << ", " << cluster.center_y() << ") ";
        os << "[";
        for (const auto& id : cluster.membersId){
            os << id << ",";
        }
        os << "]";
        return os;       
    }
};

// Wynik findClosestPair: indeksy klastrów w wektorze i kwadrat odległości między ich środkami.

struct ClosestPair{
    std::size_t i, j;
    double dist2;
};


/**
 * Parsuje linię zawierającą dokładnie jedną wartość typu T (n albo m).
 * Szablon pozwala obsłużyć liczbę całkowitą i rzeczywistą jedną funkcją.
 * @return false, gdy brakuje wartości, ma zły typ albo po niej są dodatkowe dane.
 */

template <typename T>
static bool parseSingle(const std::string& line, T& out){
    std::istringstream iss(line);
    if(!(iss >> out)){
        return false; 
    }  
    iss >> std::ws;                    
    return iss.eof();                  
}

/**
 * Parsuje linię pachołka w formacie "<id> <x> <y>".
 * @return false przy brakującym polu, złym typie lub nadmiarowych danych.
 */
static bool parseConeLine(const std::string& line, Cone& out){
    std::istringstream iss(line);
    if(!(iss >> out.id >> out.x >> out.y)){
        return false;
    }
    iss >> std::ws;
    return iss.eof();
}

/**
 * Wczytuje plik wejściowy. Struktura funkcji odpowiada strukturze pliku:
 * linia 1 - n, linia 2 - m, potem dokładnie n linii z pachołkami.
 * @throws std::runtime_error z numerem błędnej linii.
 */

InputData readInput(const std::string& filepath){
    std::ifstream file(filepath);
    if(!file) {
        throw std::runtime_error("Unable to open a file " + filepath);
    }

    std::string line;
    InputData data{};
    
    long nTemp{};
    if(!std::getline(file, line) || !parseSingle(line, nTemp) || nTemp < 0){
        throw std::runtime_error("Line 1: Wrong n number");
    }
    data.n = static_cast<std::size_t>(nTemp);

    if(!std::getline(file, line) || !parseSingle(line, data.minDist) || data.minDist <= 0){
        throw std::runtime_error("Line 2: Wrong minDist number");
    }

    data.cones.reserve(data.n);

    for(std::size_t i = 0; i < data.n; ++i){
        Cone c{};
        if(!std::getline(file, line) || !parseConeLine(line, c)){
            throw std::runtime_error("Line " + std::to_string(i+3) + ": Wrong number");
        }
        data.cones.push_back(c);
    }

    return data;
}

/**
 * Znajduje parę klastrów o najmniejszej odległości między środkami.
 * Wywoływana w każdej iteracji mergeAll, bo po scaleniu powstaje nowy środek,
 * więc odległości trzeba liczyć od nowa. 
 * Złożoność: O(k^2) dla k klastrów.
 * @return std::nullopt, gdy klastrów jest mniej niż dwa.
 */

std::optional<ClosestPair> findClosestPair(const std::vector<Cluster>& clusters){
    
    if(clusters.size() < 2){
        return std::nullopt;
    }
    // Start od infinity, żeby pierwsza sprawdzona para zawsze ją pobiła.
    ClosestPair best{0, 0, std::numeric_limits<double>::max()};
    
    for(std::size_t i = 0; i < clusters.size(); ++i){
        // j od i+1: każda para sprawdzana raz, bez porównania klastra z samym sobą
        for(std::size_t j = i+1; j < clusters.size(); ++j){
            double dx = clusters[i].center_x() - clusters[j].center_x();
            double dy = clusters[i].center_y() - clusters[j].center_y();

            double d2 = dx*dx + dy*dy;

            if(d2 < best.dist2){
                best = {i, j, d2};
            }
        }
    }

    return best;
}

/**
 * Deduplikacja pachołków.
 * Na starcie każdy pachołek jest osobnym klastrem. W każdym kroku scalana jest
 * najbliższa para klastrów, aż najbliższa para będzie w odległości >= minDist.
 * Odległość równa minDist jest dopuszczalna, więc wtedy nie scalamy.
 * Po zakończeniu żadne dwa klastry nie są bliżej niż minDist.
 * Złożoność: O(n^3).
 */

std::vector<Cluster> mergeAll(const std::vector<Cone>& cones, double minDist){
    std::vector<Cluster> clusters;
    for(const auto& cone : cones){
        clusters.push_back(Cluster{cone.x, cone.y, {cone.id}, cone.id});
    } 
    while(true){
        auto pair = findClosestPair(clusters);
        if(!pair || pair->dist2 >= minDist*minDist){
            break;
        }
        Cluster& a = clusters[pair->i];
        Cluster& b = clusters[pair->j];
        a.sumX += b.sumX;
        a.sumY += b.sumY;
        a.membersId.insert(a.membersId.end(), b.membersId.begin(), b.membersId.end());
        if (pair->j != clusters.size() - 1) {
            clusters[pair->j] = std::move(clusters.back());
        }
        clusters.pop_back();
    }

    return clusters;
}

/**
 * Iloczyn wektorowy (B - A) x (C - B), równoważny standardowemu (B - A) x (C - A),
 * bo różnią się o (B - A) x (B - A) = 0.
 * Mówi, w którą stronę skręca ścieżka A -> B -> C:
 *  > 0 - w lewo,
 *  < 0 - w prawo,
 *  = 0 - punkty współliniowe.
 */

double crossProduct(const Cone& A, const Cone& B, const Cone& C){
    return (B.x - A.x)*(C.y - B.y) - (B.y - A.y)*(C.x - B.x);
}


/**
 * Otoczka wypukła algorytmem Andrew's monotone chain, O(n log n).
 * Punkty sortujemy po x (potem y), budujemy dolną otoczkę od lewej do prawej
 * i górną od prawej do lewej, za każdym razem odrzucając punkty, w których
 * ścieżka nie skręca w lewo.
 *
 * Uwagi:
 *  - warunek <= 0 odrzuca też punkty współliniowe leżące na krawędzi otoczki
 *    (obwód się nie zmienia, nie trafiają tylko na listę ID),
 *  - po deduplikacji na wejściu nie ma zduplikowanych punktów.
 * @return wierzchołki otoczki w kolejności przeciwnej do ruchu wskazówek zegara.
 */

std::vector<Cone> getBorderPath(std::vector<Cone> cones){
    int n = cones.size();
    int k = 0;
    if (n<=3){
        return cones;
    }
    
    std::vector<Cone> hull(2*n);

    std::sort(cones.begin(), cones.end());
    for(int i = 0; i < n; ++i){
        // Usuwamy punkty, które nie skręcają w lewo
        while(k >= 2 && crossProduct(hull[k-2], hull[k-1], cones[i]) <= 0){
            k--;
        }
        hull[k++] = cones[i];
    }

    int lower_hull_size = k;

    for(int i = n-2; i >=0; --i){
        while(k >= lower_hull_size + 1 && crossProduct(hull[k-2], hull[k-1], cones[i]) <= 0){
            k--;
        }
        hull[k++] = cones[i];
    }

    hull.resize(k-1);

    return hull;
}

/// Odległość euklidesowa między dwoma pachołkami.
double calculate_distance(const Cone& a, const Cone& b){
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

/**
 * Obwód wielokąta: suma długości boków, łącznie z bokiem domykającym
 * (ostatni -> pierwszy, stąd indeks (i + 1) % n).
 */
double calculate_perimeter(const std::vector<Cone>& hull){
    double perimeter = 0.0;
    int n = hull.size();
    if (n < 2){
        return 0.0;
    }
    
    for (int i = 0; i < n; ++i) {
        perimeter += calculate_distance(hull[i], hull[(i + 1) % n]);
    }
    return perimeter;
}

/**
 * Nadaje ID klastrom po zakończeniu deduplikacji.
 *  - klaster z jednym członkiem zachowuje oryginalne ID,
 *  - scalony klaster dostaje kolejne nowe ID od firstNewId.
 * ID nadajemy dopiero po deduplikacji, żeby scalenia pośrednie nie zużywały numerów.
 * Klastry numerujemy w kolejności najmniejszego ID członka, więc wynik nie zależy
 * od kolejności linii w pliku ani od przestawień swap-and-pop.
 * @param firstNewId powinno wynosić max(id) + 1
 */

void assignIds(std::vector<Cluster>& clusters, int firstNewId){

    for(auto& c : clusters){
        std::sort(c.membersId.begin(), c.membersId.end());
    }

    std::sort(clusters.begin(), clusters.end(), 
        [](const Cluster& a, const Cluster& b){
            return a.membersId[0] < b.membersId[0];
        });

    int nextId = firstNewId;
    for(auto& c : clusters){
        if(c.membersId.size() > 1){
            c.id = nextId++;
        }
        else{
            c.id = c.membersId[0];
        }
    }
}

/// Zamienia klastry na pachołki (ID + środek) - wejście dla getBorderPath.
std::vector<Cone> toCones(const std::vector<Cluster>& clusters){
    std::vector<Cone> result;
    result.reserve(clusters.size());
    for (const auto& c : clusters) {
        result.push_back({c.id, c.center_x(), c.center_y()});
    }
    return result;
}

int main(int argc, char* argv[]){
    if(argc != 2){
        std::cerr << "Usage: " << argv[0] << " <file.txt>" << std::endl;
        return 1;
    }
    try{
        InputData input = readInput(argv[1]);
        std::cout << input << std::endl;
        std::cout << "===== MERGED =====" << std::endl;
        std::vector<Cluster> merged = mergeAll(input.cones, input.minDist);
        
        int maxId = std::max_element(input.cones.begin(), input.cones.end(),
                [](const Cone& a, const Cone& b) { return a.id < b.id; })->id;
        
        assignIds(merged, maxId + 1);

        for(const auto& cluster : merged){
            if(cluster.membersId.size() > 1){
                std::cout << cluster << std::endl;
            }
        }
        std::cout << "===== ENCLOSURE =====" << std::endl;
        std::vector<Cone> hull = getBorderPath(toCones(merged));
        std::cout << "Perimeter: " << calculate_perimeter(hull);
        std::cout << " [";
        for(const auto& c : hull){
            std::cout << c.id << ", ";
        }
        std::cout << "]" << std::endl;
    } catch(const std::exception& e){
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}


