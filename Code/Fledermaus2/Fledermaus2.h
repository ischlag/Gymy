#define _USE_MATH_DEFINES

#include <amp.h>
#include <iostream>
#include <chrono>
#include <amp_math.h>
#include <queue>
#include <random>
#include <algorithm>
#include <ctime>

#include "engine.h"

using namespace std;
using namespace concurrency;
using namespace concurrency::precise_math;

class Fledermaus2
{
public:

	//********************* GYMY BUILDERs
	void Fledermaus2::gymy_builder_serial(vector<double>& gArr, const int depth, const int size,
		vector<vector<double>>& inputArr, const int rows, const int elements, double dist[]);
	void Fledermaus2::gymy_builder_amp(vector<double>& gArr, const int depth, const int size,
		vector<vector<double>>& inputArr, const int rows, const int elements, double dist[]);

	//********************* RECONSTRUCTION 
	void Fledermaus2::rek_serial(vector<double>& gArr, const int size, const int depth,
		vector<double>& tArr, const int rows, const int elements,
		vector<double>& rArr, double dist[], int typeDist);
	void Fledermaus2::rek_amp(vector<double>& gArr, const int size, const int depth,
		vector<double>& tArr, const int rows, const int elements,
		vector<double>& rArr, double dist[], int typeDist);

	//********************* DATA PRODUCTIONS
	vector<double> Fledermaus2::get_empty_gymy(const int size, const int depth);
	vector<double> Fledermaus2::get_empty_res(const int amount, const int size);
	vector<double> Fledermaus2::get_bat_train_data();
	vector<double> Fledermaus2::get_bat_test_data();
	vector<vector<double>> Fledermaus2::get_bat_data();
	vector<vector<double>> Fledermaus2::get_iris_data();
	vector<vector<double>> Fledermaus2::generate_data(int classes, int elements, int samplesPerClass);

	//********************* UTILS
	void Fledermaus2::visualize(double* gymy, int size);
	void Fledermaus2::gymy_dig(vector<double>& gymy, int depth, int size, double threshold);
	double Fledermaus2::arr_mean_value(vector<double>& arr, const int depth, const int size);
	void Fledermaus2::print_first_x(vector<double>& arr, int x);
	int Fledermaus2::get_max_value_index(vector<double> arr, int begin, int size);
	double Fledermaus2::rek_quali(vector<double> arr, int begin, int size);

	//********************* MAIN PROGRAM
	void Fledermaus2::run();
	void Fledermaus2::run_bat();
	void Fledermaus2::run_iris();
};