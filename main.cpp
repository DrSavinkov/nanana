#include <cmath>
#include <iostream>
#include <vector>
using namespace std;
double rnd() { return rand() / (double)RAND_MAX; }
double f(double x) { return x; }
class Neyron {
  double free, free_save;
  vector<double> nonfree, nonfree_save;
  double result;
  void save() { free_save = free, nonfree_save = nonfree; }

public:
  Neyron() { free = rnd() * 2 - 1; }
  void conf(const int n) {
    for (int i = 0; i < n; i++)
      nonfree.push_back(rnd() * 2 - 1);
  }
  void change(const double d) {
    save();
    free = free + d * (rnd() * 2 - 1);
    for (auto &i : nonfree)
      i = i + d * (rnd() * 2 - 1);
  }
  void restore() { free = free_save, nonfree = nonfree_save; }
  void solve(const vector<Neyron> &prev_layer) {
    result = free;
    for (int i = 0; i < prev_layer.size(); i++)
      result += nonfree[i] * prev_layer[i].result;
    result = f(result);
  }
  double get_result() const { return result; }
  void set_result(const double v) { result = v; }
};
class Network {
  vector<vector<Neyron>> layers;

public:
  void change(const double d) {
    for (auto &i : layers)
      for (auto &j : i)
        j.change(d);
  }
  void restore() {
    for (auto &i : layers)
      for (auto &j : i)
        j.restore();
  }
  void conf(const vector<int> &size) {
    layers.resize(size.size());
    for (int i = 0; i < size.size(); i++) {
      layers[i].resize(size[i]);
      if (i)
        for (int j = 0; j < size[i]; j++)
          layers[i][j].conf(size[i - 1]);
    }
  }
  vector<double> solve(const vector<double> &input) {
    for (int i = 0; i < input.size(); i++)
      layers[0][i].set_result(input[i]);
    for (int i = 1; i < layers.size(); i++) {
      for (int j = 0; j < layers[i].size(); j++)
        layers[i][j].solve(layers[i - 1]);
    }
    vector<double> r;
    for (const auto &i : layers[layers.size() - 1])
      r.push_back(i.get_result());
    return r;
  }
};
double error(Network N, const vector<vector<double>> &input,
             const vector<vector<double>> &output) {
  double sum = 0;
  int i, j;
  for (i = 0; i < input.size(); i++) {
    auto ans = N.solve(input[i]);
    for (j = 0; j < ans.size(); j++)
      sum += pow(ans[j] - output[i][j], 2.);
  }
  return sum;
}
// double error(Network N, const vector<vector<double>> &input,
//              const vector<vector<double>> &output) {
//   double sum = 0, div, diff;
//   int i, j;
//   for (i = 0; i < input.size(); i++) {
//     auto ans = N.solve(input[i]);
//     for (j = 0; j < ans.size(); j++) {
//       diff = ans[j] - output[i][j];
//       div = fmax(fabs(ans[j]), fabs(output[i][j]));
//       if (div < 1e-16) // защита от деления на 0
//         div = 1e-16;
//       sum += pow(diff / div, 2.);
//     }
//   }
//   return sum;
// }
int main() {
  vector<vector<double>> in = {{1, 1}, {0., 1.}, {1., 0.}, {1., 1.}};
  vector<vector<double>> out = {{1}, {0}, {0.}, {0.}};
  int N = 1000;
  Network net;
  net.conf({2, 3, 1});
  double err = error(net, in, out);
  for (int i = 0; i < N; i++) {
    net.change(0.01);
    double nerr = error(net, in, out);
    if (nerr > err) {
      net.restore();
    } else {
      err = nerr;
    }
    cout << i + 1 << ") " << err << endl;
  }
  for (int i = 0; i < in.size(); i++) {
    cout << "In: (";
    for (int j = 0; j < in[i].size(); j++)
      cout << (j ? ", " : "") << in[i][j];
    cout << "): (";
    auto ans = net.solve(in[i]);
    for (int j = 0; j < ans.size(); j++)
      cout << (j ? ", " : "") << ans[j];
    cout << ") (";
    for (int j = 0; j < out[i].size(); j++)
      cout << (j ? ", " : "") << out[i][j];
    cout << ")\n";
  }
  return 0;
}