#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>

// beste Compileroptionen: -fopenmp -ffast-math -O3
// -ffast-math nicht aktivieren, wenn Sie darauf angewiesen sind,
// dass inf oder NaN explizit erkannt werden kann
// -march=native scheint hier aus unbekannten Gründen kontraproduktiv

struct Point { double x, y, z; };

int main()
{
  //
  auto StartTime = std::chrono::steady_clock::now();

  // Konstanten, Massen, Positionen initialisieren
  const size_t N = 20000;
  const double G = 6.674e-11;
  std::vector<Point> r(N);
  std::vector<double> m(N);
  for (int i = 0; i < N; ++i)
  {
    m[i] = i; // Massen willkürlich
    r[i] = {double(i), double(i), double(i)}; // Positionen willkürlich
  }

  // Matrix erstellen und berechnen
  std::vector<std::vector<double>> Matrix(N, std::vector<double>(N, 0.0));
  // parallelisiert man stattdessen die innere Schleife,
  // wird es interessanterweise noch ein kleines bisschen schneller,
  // allerdings ist Ihr N nicht so groß wie hier
  #pragma omp parallel for
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < i; ++j) {
      Matrix[i][j] = G*m[i]*m[j]/std::pow((r[i].x - r[j].x)*(r[i].x - r[j].x) + (r[i].y - r[j].y)*(r[i].y - r[j].y) + (r[i].z - r[j].z)*(r[i].z - r[j].z), 1.5);
    }
  }

  auto Duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - StartTime);
  std::cout << "Dauer: " << Duration.count() << " ms" << std::endl;

  // sinnlose Summe bilden und ausgeben, damit nichts wegoptimiert wird
  double Sum = 0.0;
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      Sum += Matrix[i][j];
    }
  }
  std::cout << "Sum = " << Sum << std::endl;

  return 0;
}
