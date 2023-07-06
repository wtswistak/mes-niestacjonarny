#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
struct Node {
    double x;
    int BC;       // 1 konwe 2 strumien ciepla
    double t;     // temp na wezle
};
struct Element {
    int id[2];       // E1-N1,N2-id[1,2]
    double Le;       // dlugosc elemn
    double** KeLocal; // lokalna macierz H
    double* FeLocal;
    Node* nodes;     // wskaznik do tablicy węzłów
};

struct GlobalData {
    double rMin;
    double rMax;
    double alfaAir;
    double tempAir;
    double tempBegin;
    double c;
    double ro;
    double k;
    double tauMax;
    int nH;
    int nE;
    int nP;
    double e1;
    double e2;
    double dTau;
    double dr;
    int w[2];

};
struct SiatkaMes {
    Element* elements;
    Node* nodes;
};

struct SOE {
    double** KeGlobal;
    double* FeGlobal;
    double* T;
};
void calcNodeElement(SiatkaMes& siatka, const GlobalData& data) {
    siatka.nodes = new Node[data.nH];
    for (int i = 0; i < data.nH; i++) {
        siatka.nodes[i].x = data.dr * i;
        siatka.nodes[i].BC = 0;
        if (i == data.nH - 1) {
            siatka.nodes[i].BC = 1;
        }
        siatka.nodes[i].t = data.tempBegin;
    }

    siatka.elements = new Element[data.nE];
    for (int i = 0; i < data.nE; i++) {
        siatka.elements[i].id[0] = i;
        siatka.elements[i].id[1] = i + 1;
        siatka.elements[i].Le = data.dr;
        siatka.elements[i].KeLocal = new double* [2];
        for (int j = 0; j < 2; j++) {
            siatka.elements[i].KeLocal[j] = new double[2];
        }

        siatka.elements[i].FeLocal = new double[2];
        siatka.elements[i].nodes = new Node[2];
        siatka.elements[i].nodes[0] = siatka.nodes[i];
        siatka.elements[i].nodes[1] = siatka.nodes[i + 1];
    }
}
void calcMacierzWektorLokalny(Element& element, const GlobalData& data) {
    double N1[2];
    double N2[2];
    N1[0] = 0.5 * (1 - data.e1);
    N1[1] = 0.5 * (1 - data.e2);
    N2[0] = 0.5 * (1 + data.e1);
    N2[1] = 0.5 * (1 + data.e2);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            element.KeLocal[i][j] = 0.0;
        }
    }

    element.FeLocal[0] = 0;
    element.FeLocal[1] = 0;

    for (int ip = 0; ip < data.nP; ip++) {
        double Rp = N1[ip] * element.nodes[0].x + N2[ip] * element.nodes[1].x;
        double TpTau = N1[ip] * element.nodes[0].t + N2[ip] * element.nodes[1].t;

        double K = data.k;
        double dR = data.dr;
        double C = data.c;
        double Ro = data.ro;
        double dTau = data.dTau;

        element.KeLocal[0][0] += K * Rp * data.w[ip] / dR + C * Ro * dR * Rp * data.w[ip] * pow(N1[ip], 2) / dTau;
        element.KeLocal[0][1] -= K * Rp * data.w[ip] / dR - C * Ro * dR * Rp * data.w[ip] * N1[ip] * N2[ip] / dTau;
        element.KeLocal[1][0] = element.KeLocal[0][1];
        element.KeLocal[1][1] += K * Rp * data.w[ip] / dR + C * Ro * dR * Rp * data.w[ip] * pow(N2[ip], 2) / dTau;

        element.FeLocal[0] -= C * Ro * dR * TpTau * Rp * data.w[ip] * N1[ip] / dTau;
        element.FeLocal[1] -= C * Ro * dR * TpTau * Rp * data.w[ip] * N2[ip] / dTau;
    }

    if (element.nodes[1].BC == 1) {
        element.FeLocal[1] -= 2 * data.alfaAir * data.rMax * data.tempAir;
        element.KeLocal[1][1] += 2 * data.alfaAir * data.rMax;
    }
}
void calcMacierzWektorGlobalny(const SiatkaMes& siatka, const GlobalData& data, SOE& soe) {
    soe.KeGlobal = new double* [data.nH];
    for (int i = 0; i < data.nH; i++) {
        soe.KeGlobal[i] = new double[data.nH];
        for (int j = 0; j < data.nH; j++) {
            soe.KeGlobal[i][j] = 0.0;
        }
    }

    soe.FeGlobal = new double[data.nH];
    for (int i = 0; i < data.nH; i++) {
        soe.FeGlobal[i] = 0.0;
    }

    for (int i = 0; i < data.nE; i++) {
        const Element& element = siatka.elements[i];
        soe.KeGlobal[element.id[0]][element.id[0]] += element.KeLocal[0][0];
        soe.KeGlobal[element.id[0]][element.id[1]] += element.KeLocal[0][1];
        soe.KeGlobal[element.id[1]][element.id[0]] += element.KeLocal[1][0];
        soe.KeGlobal[element.id[1]][element.id[1]] += element.KeLocal[1][1];

        soe.FeGlobal[element.id[0]] += element.FeLocal[0];
        soe.FeGlobal[element.id[1]] += element.FeLocal[1];
    }
}
void gauss(SOE& soe, int n) {
    double** tab = soe.KeGlobal;
    double* rozw = soe.FeGlobal;
    double* temp = new double[n];

    for (int j = 0; j < n - 1; j++) {
        for (int i = j + 1; i < n; i++) {
            double mnoznik = tab[i][j] / tab[j][j];
            for (int k = j; k < n; k++) {
                tab[i][k] = tab[i][k] - tab[j][k] * mnoznik;
            }
            rozw[i] = rozw[i] - rozw[j] * mnoznik;
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        temp[i] = rozw[i];
        for (int j = i + 1; j < n; j++) {
            temp[i] -= tab[i][j] * temp[j];
        }
        temp[i] /= tab[i][i];
    }

    for (int i = 0; i < n; i++) {
        soe.T[i] = abs(temp[i]);
    }

    delete[] temp;
}


int main() {
    GlobalData data = {
          0, // rMin
          0.08, // rMax
          300, // alfaAir
          1200, // tempAir
          100, // tempBegin
          700, // c
          7800, // ro
          25, // k
          1000, // tauMax
          8,// nH liczba wezlow
          7, // nE liczb el
          2, // nP liczba pkt calkowania
          -0.5773502692, // e1
          0.5773502692, // e2
          10 , //dtau
          (data.rMax - data.rMin) / data.nE, //dr
          { 1,1 },
    };
    SiatkaMes siatka;
    SOE soe;
    siatka.elements = new Element[data.nE];

    calcNodeElement(siatka, data);
    for (int i = 0; i < data.nE; i++) {
        calcMacierzWektorLokalny(siatka.elements[i], data);
    }
    calcMacierzWektorGlobalny(siatka, data, soe);

    cout << "Globalna Macierz Ke:" << endl;
    for (int i = 0; i < data.nH; i++) {
        for (int j = 0; j < data.nH; j++) {
            cout << soe.KeGlobal[i][j] << " ";
        }
        cout << endl;
    }

    cout << "\nGlobalny Wektor Fe:" << endl;
    for (int i = 0; i < data.nH; i++) {
        cout << soe.FeGlobal[i] << " ";
    }
    cout << endl;

    soe.T = new double[data.nH];
    gauss(soe, data.nH);

    cout << "\nWynik" << endl;
    for (int i = 0; i < data.nH; i++) {
        cout << "Temperatura " << i + 1 << ": " << soe.T[i] << endl;
    }

    delete[] siatka.nodes;
    for (int i = 0; i < data.nE; i++) {
        delete[] siatka.elements[i].KeLocal;
        delete[] siatka.elements[i].FeLocal;
        delete[] siatka.elements[i].nodes;
    }
    delete[] siatka.elements;

    return 0;
}
