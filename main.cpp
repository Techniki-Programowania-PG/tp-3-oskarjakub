#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <vector>
#include <complex>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;
using namespace matplot;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("sinus", [](double czestotliowsc, double start, double koniec, int probki, double amplituda) {
        std::vector<double> x;
        std::vector<double> y;
        x = matplot::linspace(start, koniec, probki);
        for (int i = 0; i < x.size(); i++) {
            y.push_back(std::sin(czestotliowsc * x[i]) * amplituda);
        }
        matplot::plot(x, y);
        matplot::grid(matplot::on);
        matplot::xlabel("Os X");
        matplot::ylabel("Os Y");
        matplot::show();
        matplot::save("Sinus", "png");

        return y;
        });
    m.def("cosinus", [](double czestotliowsc, double start, double koniec, int probki, double amplituda) {
        std::vector<double> x;
        std::vector<double> y;
        x = matplot::linspace(start, koniec, probki);
        for (int i = 0; i < x.size(); i++) {
            y.push_back(std::cos(czestotliowsc * x[i]) * amplituda);
        }
        matplot::plot(x, y);
        matplot::grid(matplot::on);
        matplot::xlabel("Os X");
        matplot::ylabel("Os Y");
        matplot::show();
        matplot::save("Cosinus", "png");

        return y;
        });
    m.def("prostokat", [](double czestotliowsc, double start, double koniec, int probki, double amplituda) {
        std::vector<double> x;
        std::vector<double> y;
        x = matplot::linspace(start, koniec, probki);
        for (int i = 0; i < x.size(); i++) {
            if (std::cos(czestotliowsc * x[i]) * amplituda >= 0) {
                y.push_back(amplituda);
            }
            else {
                y.push_back(-amplituda);
            }
        }
        matplot::plot(x, y);
        matplot::grid(matplot::on);
        matplot::xlabel("Os X");
        matplot::ylabel("Os Y");
        matplot::show();
        matplot::save("prostokat", "png");
        });
    m.def("pila", [](double czestotliowsc, double start, double koniec, int probki, double amplituda) {
        std::vector<double> x;
        std::vector<double> y;
        x = matplot::linspace(start, koniec, probki);
        for (int i = 0; i < x.size(); i++) {
            y.push_back((fmod((amplituda * x[i] * czestotliowsc),2))-1);
        }
        matplot::plot(x, y);
        matplot::grid(matplot::on);
        matplot::xlabel("Os X");
        matplot::ylabel("Os Y");
        matplot::show();
        matplot::save("pila", "png");
        });

    m.def("DFT", [](std::vector<double> entry, double czas_p, double czas_k, int N) {
        std::vector<double> x = matplot::linspace(czas_p, czas_k, N);
        std::vector<std::complex<double>> y;
        std::vector<double> a;
        std::complex<double> i(0.0, 1.0);
        std::complex<double> z;
        double tmp;

        for (int j = 0; j < N; j++) {
            z = 0;
            for (int k = 0; k < N; k++) {
                tmp = 2.0 * matplot::pi * double(j) * double(k);
                z += entry[k] * std::exp((-i * tmp) / double(N));
            }
            y.push_back(z);
            a.push_back(std::abs(z));
        }
        
        matplot::plot(x, a);
        matplot::grid(matplot::on);
        matplot::xlabel("Os X");
        matplot::ylabel("Os Y");
        matplot::show();
        matplot::save("dft", "png");

        return y;

        });
    m.def("IDFT", [](std::vector<std::complex<double>> entry, double czas_p, double czas_k, int N) {
        std::vector<double> x = matplot::linspace(czas_p, czas_k, N);
        std::vector<std::complex<double>> y;
        std::vector<double> a;
        std::complex<double> i(0.0, 1.0);
        std::complex<double> z;
        double tmp;

        for (int j = 0; j < N; j++) {
            z = 0;
            for (int k = 0; k < N; k++) {
                tmp = 2.0 * matplot::pi * double(j) * double(k);
                z += entry[k] * std::exp((i * tmp) / double(N));
            }
            z /= double(N);
            y.push_back(z);
            a.push_back(std::real(z));
        }

        matplot::plot(x, a);
        matplot::grid(matplot::on);
        matplot::xlabel("Os X");
        matplot::ylabel("Os Y");
        matplot::show();
        matplot::save("idft", "png");

        return y;
        });

    m.def("unoD", [](std::vector<double> entry, double czas_p, double czas_k, int N) {
        std::vector<double> y;
        std::vector<double> x = matplot::linspace(czas_p, czas_k, N);
        std::vector<double> f;

        double a = 2;
        double sum = 0.0;

        for (int i = 0; i < entry.size(); i++) {
            f.push_back((pow(-1, i + 1)) * a);
            a *= 2;
        }
        for (int i = 0; i < entry.size(); i++) {
            sum = 0.0;

            for (int j = 0; j < i + 1; j++) {
                sum += entry[j] * f[i - j];
            }
            y.push_back(sum);
        }

        matplot::plot(x, y);
        matplot::xlabel("Os x");
        matplot::ylabel("Os y");
        matplot::grid(matplot::on);
        matplot::show();
        matplot::save("1d", "png");
    });

    m.def("dosD", [](double czas_p, double czas_k) {
        czas_k += 1;
        int size = czas_k + 1;
        double wartosc = 1.0;

        std::vector<double> x = linspace(czas_p, czas_k, size);
        std::vector<double> y = linspace(czas_p, czas_k, size);
        std::vector<std::vector<double>> X(size, std::vector<double>(size));
        std::vector<std::vector<double>> Y(size, std::vector<double>(size));
        std::vector<std::vector<double>> dane(size, std::vector<double>(size));
        std::vector<std::vector<double>> f{
            {0, -1, 0},
            {-1, 5, -1},
            {0, -1, 0}
        };
        std::vector<std::vector<double>> Z(size, std::vector<double>(size, 0.0));

        for (int i = 0; i < size - 1; ++i) {
            for (int j = 0; j < size - 1; ++j) {
                X[i][j] = x[j];
                Y[i][j] = y[i];
            }
        }
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                dane[i][j] = wartosc;
                wartosc += 0.5;
            }
        }
        for (int i = 1; i < size - 1; ++i) {
            for (int j = 1; j < size - 1; ++j) {
                double suma = 0.0;
                for (int k = -1; k <= 1; ++k) {
                    for (int m = -1; m <= 1; ++m) {
                        suma += f[k + 1][m + 1] * dane[i + k][j + m];
                    }
                }
                Z[i][j] = suma;
            }
        }

        surf(X, Y, Z);
        xlabel("Os x");
        ylabel("Os y");
        zlabel("Os z");
        grid(on);
        show();
        matplot::save("2d", "png");

        return Z;
        });

    m.def("Piki", [](std::vector<double> entry, double czas_p, double czas_k, int N) {
        std::vector<double> piki_nr;
        std::vector<double> piki;
        std::vector<double> x = matplot::linspace(czas_p, czas_k, N);
        
        for (int i = 1; i < entry.size() - 1; i++) {
            if (entry[i] > entry[i - 1] && entry[i] > entry[i + 1])
            {
                piki_nr.push_back(i);
            }
        }
        for (int i = 0; i < piki_nr.size(); i++)
        {
            piki.push_back(x[piki_nr[i]]);
        }

        return piki;
        });

    m.attr("__version__") = "dev";
}
