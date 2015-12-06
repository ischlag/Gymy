#include "Fledermaus2.h"

const int step = 4;

//********************* GYMY BUILDERs
void Fledermaus2::gymy_builder_serial(vector<double>& gArr, const int depth, const int size,
	vector<vector<double>>& inputArr, const int rows, const int elements, double dist[]) {
	const int q = 1;

	for (int x2 = 0; x2 < depth; x2++) {
		for (int x1 = 0; x1 < size; x1++) {
			double sum1 = 0;
			for (int j = 0; j < rows; j++) {
				double sum2 = 0;
				for (int i = 0; i < elements; i++) {
					double tmp = fabs(cos(2 * M_PI / q * sqrt(pow(dist[i], 2) + pow(x1*step - inputArr[j][i], 2) + pow(x2*step, 2))));
					sum2 += tmp < 0.707 ? -1 : 1;
					//sum2 += (tmp < 0.95 ? -0.25 : 1);
				}
				sum1 += sum2*sum2;
			}
			gArr[x2 * size + x1] = sum1;
		}
	}
}
void Fledermaus2::gymy_builder_amp(vector<double>& gArr, const int depth, const int size,
	vector<vector<double>>& inputArr, const int rows, const int elements, double dist[]) {
	const int q = 1;
	
	vector<double>* currData = new vector<double>();
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < elements; j++) {
			currData->push_back(inputArr[i][j]);
		}
	}

	array_view<double, 2> gA(depth, size, gArr);
	array_view<const double, 2> input(rows, elements, *currData);
	array_view<const double> d(elements, dist);
		
	//for (int x2 = 0; x2 < depth; x2++) {
		//for (int x1 = 0; x1 < size; x1++) {
	parallel_for_each(
		gA.extent,
		[=](index<2> idx) restrict(amp)
		{		
			double x2 = idx[0], x1 = idx[1];
			double sum1 = 0;
			for (int j = 0; j < rows; j++) {
				double sum2 = 0;
				for (int i = 0; i < elements; i++) {
					double tmp = fabs(cos(2 * M_PI / q * sqrt(pow(d[i], 2) + pow(x1*step - input(j,i), 2) + pow(x2*step, 2))));
					sum2 += tmp < 0.707 ? -1 : 1;
					//sum2 += (tmp < 0.95 ? -0.25 : 1);
				}
				sum1 += sum2*sum2;
			}
			gA[idx] = sum1;			
		});
	gA.synchronize();
		//}
	//}
}

//********************* RECONSTRUCTION 
void Fledermaus2::rek_serial(vector<double>& gArr, const int size, const int depth,
	vector<double>& tArr, const int rows, const int elements,
	vector<double>& rArr, double dist[], int typeDist){
	const int q = 1;

	for (int y = 0; y < size; y++) {
		double sum = 0;
		for (int i = 0; i < elements; i++) {
			for (int x2 = 0; x2 < depth; x2++) {
				for (int x1 = 0; x1 < size; x1++) {
					double tmp = fabs(cos(2 * M_PI / q*(sqrt(pow(dist[i], 2) + pow(x1*step - tArr[i], 2) + pow(x2*step, 2)) + sqrt(pow(dist[typeDist], 2) + pow(x1*step - y*step, 2) + pow(x2*step, 2)))));
					sum += gArr[x2* size + x1] * (tmp < 0.707 ? -1 : 1);
					//sum += gArr[x2* size + x1] * (tmp < 0.95 ? -0.25 : 1);
				}
			}
		}
		rArr[y] = sum*sum;
	}
}
void Fledermaus2::rek_amp(vector<double>& gArr, const int size, const int depth,
	vector<double>& tArr, const int rows, const int elements,
	vector<double>& rArr, double dist[], int typeDist){
	const int q = 1;
		
	array_view<const double, 2> gA(depth, size, gArr);
	array_view<const double> tA(elements, tArr);
	array_view<double> rA(size, rArr);
	array_view<const double> d(elements, dist);

	//for (int y = 0; y < size; y++) {
	parallel_for_each(
		rA.extent,
		[=](index<1> idx) restrict(amp)
	{
		int y = idx[0];
		double sum = 0;
		for (int i = 0; i < elements; i++) {
			for (int x2 = 0; x2 < depth; x2++) {
				for (int x1 = 0; x1 < size; x1++) {
					double tmp = fabs(cos(2 * M_PI / q*(sqrt(pow(d[i], 2) + pow(x1*step - tA(i), 2) + pow((float)x2*step, 2)) + sqrt(pow(d[typeDist], 2) + pow((float)x1*step - y*step, 2) + pow((float)x2*step, 2)))));
					sum += gA(x2,x1) * (tmp < 0.707 ? -1 : 1);
					//sum += gA(x2,x1) * (tmp < 0.95 ? -0.25 : 1);
				}
			}
		}
		rA(y) = sum*sum;
	});
	rA.synchronize();
	//}
}

//********************* DATA PRODUCTIONS
vector<double> Fledermaus2::get_empty_gymy(const int size, const int depth) {
	vector<double>* new_gymy = new vector<double>(size*depth, 0.0);
	return *new_gymy;
}
vector<double> Fledermaus2::get_empty_res(const int amount, const int size) {
	vector<double>* new_gymy = new vector<double>(size*amount, 0.0);
	return *new_gymy;
}
vector<vector<double>> Fledermaus2::get_bat_data() {
	vector<vector<double>> data{

	vector<double>{ 340, 204, 167, 264, 100},
	vector<double>{ 212, 180, 164, 371, 100},
	vector<double>{ 433, 214, 185, 226, 100},
	vector<double>{ 431, 214, 183, 235, 100},
	vector<double>{ 340, 189, 157, 286, 100},
	vector<double>{ 354, 198, 154, 294, 100},
	vector<double>{ 369, 192, 155, 298, 100},
	vector<double>{ 357, 180, 150, 307, 100},

	vector<double>{ 248, 198, 174, 386, 200},
	vector<double>{ 297, 198, 174, 338, 200},
	vector<double>{ 262, 198, 176, 322, 200},
	vector<double>{ 245, 192, 171, 354, 200},
	vector<double>{ 242, 192, 176, 354, 200},
	vector<double>{ 257, 201, 174, 327, 200},
	vector<double>{ 254, 192, 169, 399, 200},
	vector<double>{ 404, 229, 173, 369, 200},

	vector<double>{ 494, 241, 217, 329, 300},
	vector<double>{ 393, 241, 219, 247, 300},
	vector<double>{ 461, 244, 216, 344, 300},
	vector<double>{ 509, 232, 216, 271, 300},
	vector<double>{ 392, 226, 212, 282, 300},
	vector<double>{ 366, 226, 211, 288, 300},
	vector<double>{ 469, 232, 204, 275, 300},
	vector<double>{ 402, 226, 207, 270, 300},

	vector<double>{ 390, 293, 274, 258, 400},
	vector<double>{ 376, 287, 271, 252, 400},
	vector<double>{ 475, 296, 274, 249, 400},
	vector<double>{ 376, 290, 269, 363, 400},
	vector<double>{ 378, 296, 266, 378, 400},
	vector<double>{ 409, 302, 256, 361, 400},
	vector<double>{ 481, 287, 264, 385, 400},
	vector<double>{ 497, 290, 269, 359, 400},

	//Test
	vector<double>{ 343, 180, 154, 326, 100},
	vector<double>{ 390, 238, 171, 327, 200},
	vector<double>{ 428, 235, 216, 310, 300},
	vector<double>{ 409, 290, 267, 346, 400}
};
	return data;
}
vector<vector<double>> Fledermaus2::generate_data(int classes, int samplesPerClass, int elements) {
	vector<vector<double>> data;

	for (int c = 1; c <= classes; c++) {
		int currClass = c * 100;
		for (int s = 0; s < samplesPerClass; s++) {
			vector<double> sample;
			for (int e = 0; e < elements; e++) {
				double currElement = rand() % 500;
				sample.push_back(currElement);
			}
			sample.push_back(currClass);
			data.push_back(sample);
		}
	}

	return data;
}
vector<double> Fledermaus2::get_bat_train_data() {
	vector<double>* data = new vector<double>{
		340, 204, 167, 264, 100,
		212, 180, 164, 371, 100,
		433, 214, 185, 226, 100,
		431, 214, 183, 235, 100,
		340, 189, 157, 286, 100,
		354, 198, 154, 294, 100,
		369, 192, 155, 298, 100,
		357, 180, 150, 307, 100,

		248, 198, 174, 386, 200,
		297, 198, 174, 338, 200,
		262, 198, 176, 322, 200,
		245, 192, 171, 354, 200,
		242, 192, 176, 354, 200,
		257, 201, 174, 327, 200,
		254, 192, 169, 399, 200,
		404, 229, 173, 369, 200,

		494, 241, 217, 329, 300,
		393, 241, 219, 247, 300,
		461, 244, 216, 344, 300,
		509, 232, 216, 271, 300,
		392, 226, 212, 282, 300,
		366, 226, 211, 288, 300,
		469, 232, 204, 275, 300,
		402, 226, 207, 270, 300,

		390, 293, 274, 258, 400,
		376, 287, 271, 252, 400,
		475, 296, 274, 249, 400,
		376, 290, 269, 363, 400,
		378, 296, 266, 378, 400,
		409, 302, 256, 361, 400,
		481, 287, 264, 385, 400,
		497, 290, 269, 359, 400
	};
	return *data;
}
vector<double> Fledermaus2::get_bat_test_data() {
	vector<double>* data = new vector<double>{
		343, 180, 154, 326, 100,
		390, 238, 171, 327, 200,
		428, 235, 216, 310, 300,
		409, 290, 267, 346, 400
	};
	return *data;
}
vector<vector<double>> Fledermaus2::get_iris_data(){
	vector<vector<double>> data{
		vector<double>{ 5.1000, 3.5000, 1.4000, 0.2000, 100},
		vector<double>{ 4.9000, 3.0000, 1.4000, 0.2000, 100},
		vector<double>{ 4.7000, 3.2000, 1.3000, 0.2000, 100},
		vector<double>{ 4.6000, 3.1000, 1.5000, 0.2000, 100},
		vector<double>{ 5.0000, 3.6000, 1.4000, 0.2000, 100},
		vector<double>{ 5.4000, 3.9000, 1.7000, 0.4000, 100},
		vector<double>{ 4.6000, 3.4000, 1.4000, 0.3000, 100},
		vector<double>{ 5.0000, 3.4000, 1.5000, 0.2000, 100},
		vector<double>{ 4.4000, 2.9000, 1.4000, 0.2000, 100},
		vector<double>{ 4.9000, 3.1000, 1.5000, 0.1000, 100},
		vector<double>{ 5.4000, 3.7000, 1.5000, 0.2000, 100},
		vector<double>{ 4.8000, 3.4000, 1.6000, 0.2000, 100},
		vector<double>{ 4.8000, 3.0000, 1.4000, 0.1000, 100},
		vector<double>{ 4.3000, 3.0000, 1.1000, 0.1000, 100},
		vector<double>{ 5.8000, 4.0000, 1.2000, 0.2000, 100},
		vector<double>{ 5.7000, 4.4000, 1.5000, 0.4000, 100},
		vector<double>{ 5.4000, 3.9000, 1.3000, 0.4000, 100},
		vector<double>{ 5.1000, 3.5000, 1.4000, 0.3000, 100},
		vector<double>{ 5.7000, 3.8000, 1.7000, 0.3000, 100},
		vector<double>{ 5.1000, 3.8000, 1.5000, 0.3000, 100},
		vector<double>{ 5.4000, 3.4000, 1.7000, 0.2000, 100},
		vector<double>{ 5.1000, 3.7000, 1.5000, 0.4000, 100},
		vector<double>{ 4.6000, 3.6000, 1.0000, 0.2000, 100},
		vector<double>{ 5.1000, 3.3000, 1.7000, 0.5000, 100},
		vector<double>{ 4.8000, 3.4000, 1.9000, 0.2000, 100},
		vector<double>{ 5.0000, 3.0000, 1.6000, 0.2000, 100},
		vector<double>{ 5.0000, 3.4000, 1.6000, 0.4000, 100},
		vector<double>{ 5.2000, 3.5000, 1.5000, 0.2000, 100},
		vector<double>{ 5.2000, 3.4000, 1.4000, 0.2000, 100},
		vector<double>{ 4.7000, 3.2000, 1.6000, 0.2000, 100},
		vector<double>{ 4.8000, 3.1000, 1.6000, 0.2000, 100},
		vector<double>{ 5.4000, 3.4000, 1.5000, 0.4000, 100},
		vector<double>{ 5.2000, 4.1000, 1.5000, 0.1000, 100},
		vector<double>{ 5.5000, 4.2000, 1.4000, 0.2000, 100},
		vector<double>{ 4.9000, 3.1000, 1.5000, 0.1000, 100},
		vector<double>{ 5.0000, 3.2000, 1.2000, 0.2000, 100},
		vector<double>{ 5.5000, 3.5000, 1.3000, 0.2000, 100},
		vector<double>{ 4.9000, 3.1000, 1.5000, 0.1000, 100},
		vector<double>{ 4.4000, 3.0000, 1.3000, 0.2000, 100},
		vector<double>{ 5.1000, 3.4000, 1.5000, 0.2000, 100},
		vector<double>{ 5.0000, 3.5000, 1.3000, 0.3000, 100},
		vector<double>{ 4.5000, 2.3000, 1.3000, 0.3000, 100},
		vector<double>{ 4.4000, 3.2000, 1.3000, 0.2000, 100},
		vector<double>{ 5.0000, 3.5000, 1.6000, 0.6000, 100},
		vector<double>{ 5.1000, 3.8000, 1.9000, 0.4000, 100},
		vector<double>{ 7.0000, 3.2000, 4.7000, 1.4000, 200},
		vector<double>{ 6.4000, 3.2000, 4.5000, 1.5000, 200},
		vector<double>{ 6.9000, 3.1000, 4.9000, 1.5000, 200},
		vector<double>{ 5.5000, 2.3000, 4.0000, 1.3000, 200},
		vector<double>{ 6.5000, 2.8000, 4.6000, 1.5000, 200},
		vector<double>{ 5.7000, 2.8000, 4.5000, 1.3000, 200},
		vector<double>{ 6.3000, 3.3000, 4.7000, 1.6000, 200},
		vector<double>{ 4.9000, 2.4000, 3.3000, 1.0000, 200},
		vector<double>{ 6.6000, 2.9000, 4.6000, 1.3000, 200},
		vector<double>{ 5.2000, 2.7000, 3.9000, 1.4000, 200},
		vector<double>{ 5.0000, 2.0000, 3.5000, 1.0000, 200},
		vector<double>{ 5.9000, 3.0000, 4.2000, 1.5000, 200},
		vector<double>{ 6.0000, 2.2000, 4.0000, 1.0000, 200},
		vector<double>{ 6.1000, 2.9000, 4.7000, 1.4000, 200},
		vector<double>{ 5.6000, 2.9000, 3.6000, 1.3000, 200},
		vector<double>{ 6.7000, 3.1000, 4.4000, 1.4000, 200},
		vector<double>{ 5.6000, 3.0000, 4.5000, 1.5000, 200},
		vector<double>{ 5.8000, 2.7000, 4.1000, 1.0000, 200},
		vector<double>{ 6.2000, 2.2000, 4.5000, 1.5000, 200},
		vector<double>{ 5.6000, 2.5000, 3.9000, 1.1000, 200},
		vector<double>{ 5.9000, 3.2000, 4.8000, 1.8000, 200},
		vector<double>{ 6.1000, 2.8000, 4.0000, 1.3000, 200},
		vector<double>{ 6.3000, 2.5000, 4.9000, 1.5000, 200},
		vector<double>{ 6.1000, 2.8000, 4.7000, 1.2000, 200},
		vector<double>{ 6.4000, 2.9000, 4.3000, 1.3000, 200},
		vector<double>{ 6.6000, 3.0000, 4.4000, 1.4000, 200},
		vector<double>{ 6.8000, 2.8000, 4.8000, 1.4000, 200},
		vector<double>{ 6.7000, 3.0000, 5.0000, 1.7000, 200},
		vector<double>{ 6.0000, 2.9000, 4.5000, 1.5000, 200},
		vector<double>{ 5.7000, 2.6000, 3.5000, 1.0000, 200},
		vector<double>{ 5.5000, 2.4000, 3.8000, 1.1000, 200},
		vector<double>{ 5.5000, 2.4000, 3.7000, 1.0000, 200},
		vector<double>{ 5.8000, 2.7000, 3.9000, 1.2000, 200},
		vector<double>{ 6.0000, 2.7000, 5.1000, 1.6000, 200},
		vector<double>{ 5.4000, 3.0000, 4.5000, 1.5000, 200},
		vector<double>{ 6.0000, 3.4000, 4.5000, 1.6000, 200},
		vector<double>{ 6.7000, 3.1000, 4.7000, 1.5000, 200},
		vector<double>{ 6.3000, 2.3000, 4.4000, 1.3000, 200},
		vector<double>{ 5.6000, 3.0000, 4.1000, 1.3000, 200},
		vector<double>{ 5.5000, 2.5000, 4.0000, 1.3000, 200},
		vector<double>{ 5.5000, 2.6000, 4.4000, 1.2000, 200},
		vector<double>{ 6.1000, 3.0000, 4.6000, 1.4000, 200},
		vector<double>{ 5.8000, 2.6000, 4.0000, 1.2000, 200},
		vector<double>{ 5.0000, 2.3000, 3.3000, 1.0000, 200},
		vector<double>{ 5.6000, 2.7000, 4.2000, 1.3000, 200},
		vector<double>{ 6.3000, 3.3000, 6.0000, 2.5000, 300},
		vector<double>{ 5.8000, 2.7000, 5.1000, 1.9000, 300},
		vector<double>{ 7.1000, 3.0000, 5.9000, 2.1000, 300},
		vector<double>{ 6.3000, 2.9000, 5.6000, 1.8000, 300},
		vector<double>{ 6.5000, 3.0000, 5.8000, 2.2000, 300},
		vector<double>{ 7.6000, 3.0000, 6.6000, 2.1000, 300},
		vector<double>{ 4.9000, 2.5000, 4.5000, 1.7000, 300},
		vector<double>{ 7.3000, 2.9000, 6.3000, 1.8000, 300},
		vector<double>{ 6.7000, 2.5000, 5.8000, 1.8000, 300},
		vector<double>{ 7.2000, 3.6000, 6.1000, 2.5000, 300},
		vector<double>{ 6.5000, 3.2000, 5.1000, 2.0000, 300},
		vector<double>{ 6.4000, 2.7000, 5.3000, 1.9000, 300},
		vector<double>{ 6.8000, 3.0000, 5.5000, 2.1000, 300},
		vector<double>{ 5.7000, 2.5000, 5.0000, 2.0000, 300},
		vector<double>{ 5.8000, 2.8000, 5.1000, 2.4000, 300},
		vector<double>{ 6.4000, 3.2000, 5.3000, 2.3000, 300},
		vector<double>{ 6.5000, 3.0000, 5.5000, 1.8000, 300},
		vector<double>{ 7.7000, 3.8000, 6.7000, 2.2000, 300},
		vector<double>{ 7.7000, 2.6000, 6.9000, 2.3000, 300},
		vector<double>{ 6.0000, 2.2000, 5.0000, 1.5000, 300},
		vector<double>{ 6.9000, 3.2000, 5.7000, 2.3000, 300},
		vector<double>{ 5.6000, 2.8000, 4.9000, 2.0000, 300},
		vector<double>{ 7.7000, 2.8000, 6.7000, 2.0000, 300},
		vector<double>{ 6.3000, 2.7000, 4.9000, 1.8000, 300},
		vector<double>{ 6.7000, 3.3000, 5.7000, 2.1000, 300},
		vector<double>{ 7.2000, 3.2000, 6.0000, 1.8000, 300},
		vector<double>{ 6.2000, 2.8000, 4.8000, 1.8000, 300},
		vector<double>{ 6.1000, 3.0000, 4.9000, 1.8000, 300},
		vector<double>{ 6.4000, 2.8000, 5.6000, 2.1000, 300},
		vector<double>{ 7.2000, 3.0000, 5.8000, 1.6000, 300},
		vector<double>{ 7.4000, 2.8000, 6.1000, 1.9000, 300},
		vector<double>{ 7.9000, 3.8000, 6.4000, 2.0000, 300},
		vector<double>{ 6.4000, 2.8000, 5.6000, 2.2000, 300},
		vector<double>{ 6.3000, 2.8000, 5.1000, 1.5000, 300},
		vector<double>{ 6.1000, 2.6000, 5.6000, 1.4000, 300},
		vector<double>{ 7.7000, 3.0000, 6.1000, 2.3000, 300},
		vector<double>{ 6.3000, 3.4000, 5.6000, 2.4000, 300},
		vector<double>{ 6.4000, 3.1000, 5.5000, 1.8000, 300},
		vector<double>{ 6.0000, 3.0000, 4.8000, 1.8000, 300},
		vector<double>{ 6.9000, 3.1000, 5.4000, 2.1000, 300},
		vector<double>{ 6.7000, 3.1000, 5.6000, 2.4000, 300},
		vector<double>{ 6.9000, 3.1000, 5.1000, 2.3000, 300},
		vector<double>{ 5.8000, 2.7000, 5.1000, 1.9000, 300},
		vector<double>{ 6.8000, 3.2000, 5.9000, 2.3000, 300},
		vector<double>{ 6.7000, 3.3000, 5.7000, 2.5000, 300},


		//test
		vector<double>{ 4.8000, 3.0000, 1.4000, 0.3000, 100},
		vector<double>{ 5.1000, 3.8000, 1.6000, 0.2000, 100},
		vector<double>{ 4.6000, 3.2000, 1.4000, 0.2000, 100},
		vector<double>{ 5.3000, 3.7000, 1.5000, 0.2000, 100},
		vector<double>{ 5.0000, 3.3000, 1.4000, 0.2000, 100}, 
		vector<double>{ 5.7000, 3.0000, 4.2000, 1.2000, 200},
		vector<double>{ 5.7000, 2.9000, 4.2000, 1.3000, 200},
		vector<double>{ 6.2000, 2.9000, 4.3000, 1.3000, 200},
		vector<double>{ 5.1000, 2.5000, 3.0000, 1.1000, 200},
		vector<double>{ 5.7000, 2.8000, 4.1000, 1.3000, 200},
		vector<double>{ 6.7000, 3.0000, 5.2000, 2.3000, 300},
		vector<double>{ 6.3000, 2.5000, 5.0000, 1.9000, 300},
		vector<double>{ 6.5000, 3.0000, 5.2000, 2.0000, 300},
		vector<double>{ 6.2000, 3.4000, 5.4000, 2.3000, 300},
		vector<double>{ 5.9000, 3.0000, 5.1000, 1.8000, 300}
	};
	return data;
}

//********************* UTILS
void Fledermaus2::visualize(double* gymy, int size){
	Engine* ep = engOpen(NULL);

	mxArray* G = mxCreateDoubleMatrix(1, size, mxREAL);
	memcpy((void *)mxGetPr(G), (void *)gymy, size*sizeof(double));
	engPutVariable(ep, "G", G);

	engEvalString(ep, "figure('Position',[400,300,1000,300]);");
	engEvalString(ep, "bar(G, 'histc');");
	//engEvalString(ep, "set(gca,'xticklabel',{[]})");

}
void Fledermaus2::gymy_dig(vector<double>& gymy, int depth, int size, double threshold){
	//const int step = 4;
	for (int d = 0; d < depth; d++) {
		for (int s = 0; s < size; s++) {
			if (gymy[d*size + s] < threshold) {
				gymy[d*size + s] = 0;
			}
			else {
				gymy[d*size + s] = 1;
			}
		}
	}
}
double Fledermaus2::arr_mean_value(vector<double>& arr, const int depth, const int size) {
	double sum = 0;
	int counter = 0;
	for (int d = 0; d < depth; d++) {
		for (int s = 0; s < size; s++) {
			counter++;
			sum += arr[d*size + s];
		}
	}
	cout << "Gymypunkte: " << counter << endl;
	return sum / counter;
}
void Fledermaus2::print_first_x(vector<double>& arr, int x) {
	for (int i = 0; i < x; i++) {
		cout << arr[i] << "\t";
	}cout << endl;
}
double Fledermaus2::rek_quali(vector<double> arr, int begin, int size) {
	int bestIndex = -1;
	double bestValue = -1;
	int sndBestIndex = -1;
	double sndBestValue = -1;
	for (int i = begin, j = 0; i < begin + size; i++, j++) {
		if (arr[i] > bestValue) {
			sndBestValue = bestValue;
			sndBestIndex = bestIndex;
			bestValue = arr[i];
			bestIndex = j;
		}
		else if (arr[i] > sndBestValue) {
			sndBestValue = arr[i];
			sndBestIndex = j;
		}
	}
	return (bestValue-sndBestValue)/bestValue;
}
int Fledermaus2::get_max_value_index(vector<double> arr, int begin, int size) {
	int bestIndex = -1;
	double bestValue = -1;
	for (int i = begin, j = 0; i < begin + size; i++, j++) {
		if (arr[i] > bestValue) {
			bestValue = arr[i];
			bestIndex = j;
		}
	}
	return bestIndex;
}

//********************* MAIN PROGRAM
void Fledermaus2::run() {

	const int size = 250;
	const int depth = 250;

	const int classes = 10;
	const int samplesPerClass = 3;
	const int elements = 50;
	const int elem_to_ignore = 0;

	double distances[elements + 1];
	for (int i = 0; i < elements; i++) {
		distances[i] = 300 + i * 50;
	}

	clock_t begin = clock();

	//Input
	vector<vector<double>> iArr = generate_data(classes, samplesPerClass, elements);

	//Gymy bauen
	const int trainCount = classes*samplesPerClass; //Anzahl der Fledermausdaten die genutzt werden sollen. Maximal 32
	auto gArr = get_empty_gymy(size, depth);

	gymy_builder_amp   (gArr, depth, size, iArr, trainCount, elements + 1, distances);

	double mean = arr_mean_value(gArr, depth, size);

	cout << "Mean: " << mean << endl;

	print_first_x(gArr, 20);
	gymy_dig(gArr, depth, size, mean);
	print_first_x(gArr, 20);


	//Extraktion
	vector<vector<double>> data = iArr;
	
	const int dataRows = classes*samplesPerClass;
	const int dataSize = elements - elem_to_ignore;
	const int resID = elements;

	double sumQuali = 0;
	double errors = 0;

	for (int i = 0; i < dataRows; i++) {
		vector<double> tArr(dataSize);
		for (int j = 0; j < dataSize; j++)
			tArr[j] = iArr[i][j];

		cout << "\nQuery:\t";
		for (int j = 0; j < tArr.size(); j++)
			cout << tArr[j] << "\t";
		cout << endl;

		auto rArr = get_empty_res(1, size);
		rek_amp(gArr, size, depth, tArr, dataRows, dataSize, rArr, distances, resID);
		double q = rek_quali(rArr, 0, size);
		sumQuali += q;
		double res = get_max_value_index(rArr, 0, size)*step;
		cout << "Result " << i << ": " << res << "/" << iArr[i][elements] << " Quali: " << q << endl;
		if (res != iArr[i][elements])
			errors++;
		//visualize(&rArr[0], size);
	}

	sumQuali /= dataRows;
	clock_t end = clock();

	cout << "GymyPoints:\t" << size << ", " << depth << endl;
	cout << "classes:\t" << classes << endl;
	cout << "samples:\t" << classes*samplesPerClass << endl;
	cout << "elements:\t" << elements << endl;
	cout << "elementsToIgnore:\t" << elem_to_ignore << endl;
	cout << "Average Rek Quali:\t" << sumQuali << endl;
	cout << "Rek Errors:\t" << errors << endl;
	cout << "Time:\t" << double(end - begin) / CLOCKS_PER_SEC << endl;
	cout << "Fertig" << endl;

}
void Fledermaus2::run_bat() {

	const int size = 250;
	const int depth = 250;

	const int classes = 4;
	const int samplesPerClass = 8;
	const int elements = 4;
	const int elem_to_ignore = 0;

	double distances[elements + 1];
	for (int i = 0; i < elements; i++) {
		distances[i] = 300 + i * 50;
	}

	clock_t begin = clock();

	//Input
	vector<vector<double>> iArr = get_bat_data();

	//Gymy bauen
	const int trainCount = classes*samplesPerClass; //Anzahl der Fledermausdaten die genutzt werden sollen. Maximal 32
	auto gArr = get_empty_gymy(size, depth);

	//gymy_builder_serial(gArr, depth, size, iArr, trainCount, elements + 1, distances);
	gymy_builder_amp(gArr, depth, size, iArr, trainCount, elements + 1, distances);

	double mean = arr_mean_value(gArr, depth, size);

	cout << "Mean: " << mean << endl;

	print_first_x(gArr, 20);
	gymy_dig(gArr, depth, size, mean);
	print_first_x(gArr, 20);


	//Extraktion
	vector<vector<double>> data = iArr;

	const int dataRows = classes*samplesPerClass;
	const int dataSize = elements - elem_to_ignore;
	const int resID = elements;

	double sumQuali = 0;
	double errors = 0;

	for (int i = 0; i < dataRows; i++) {
		vector<double> tArr(dataSize);
		for (int j = 0; j < dataSize; j++)
			tArr[j] = iArr[i][j];

		cout << "\nQuery:\t";
		for (int j = 0; j < tArr.size(); j++)
			cout << tArr[j] << "\t";
		cout << endl;

		auto rArr = get_empty_res(1, size);
		//rek_serial(gArr, size, depth, tArr, dataRows, dataSize, rArr, distances, resID);
		rek_amp(gArr, size, depth, tArr, dataRows, dataSize, rArr, distances, resID);
		double q = rek_quali(rArr, 0, size);
		sumQuali += q;
		double res = get_max_value_index(rArr, 0, size)*step;
		cout << "Result " << i << ": " << res << "/" << iArr[i][elements] << " Quali: " << q << endl;
		if (res != iArr[i][elements])
			errors++;
		//visualize(&rArr[0], size);
	}

	sumQuali /= dataRows;
	clock_t end = clock();

	cout << "GymyPoints:\t" << size << ", " << depth << endl;
	cout << "classes:\t" << classes << endl;
	cout << "samples:\t" << classes*samplesPerClass << endl;
	cout << "elements:\t" << elements << endl;
	cout << "elementsToIgnore:\t" << elem_to_ignore << endl;
	cout << "Average Rek Quali:\t" << sumQuali << endl;
	cout << "Rek Errors:\t" << errors << endl;
	cout << "Time:\t" << double(end - begin) / CLOCKS_PER_SEC << endl;
	cout << "Fertig" << endl;

}
void Fledermaus2::run_iris() {

	const int size = 201;
	const int depth = 201;

	const int classes = 3;
	const int samplesPerClass = 45;
	const int elements = 4;
	const int elem_to_ignore = 0;

	double distances[elements + 1];
	for (int i = 0; i < elements; i++) {
		distances[i] = 300 + i * 50;
	}

	clock_t begin = clock();

	//Input
	vector<vector<double>> iArr = get_iris_data();
	for (int i = 0; i < 150; i++) {
		for (int j = 0; j < 4; j++) {
			iArr[i][j] *= 10;
		}
	}

	//Gymy bauen
	const int trainCount = classes*samplesPerClass; //Anzahl der Fledermausdaten die genutzt werden sollen. Maximal 32
	auto gArr = get_empty_gymy(size, depth);

	//gymy_builder_serial(gArr, depth, size, iArr, trainCount, elements + 1, distances);
	gymy_builder_amp(gArr, depth, size, iArr, trainCount, elements + 1, distances);

	double mean = arr_mean_value(gArr, depth, size);

	cout << "Mean: " << mean << endl;

	print_first_x(gArr, 20);
	gymy_dig(gArr, depth, size, mean);
	print_first_x(gArr, 20);


	//Extraktion
	vector<vector<double>> data = iArr;

	const int dataRows = classes*samplesPerClass;
	const int dataSize = elements - elem_to_ignore;
	const int resID = elements;

	double sumQuali = 0;
	double errors = 0;

	for (int i = 0; i < dataRows; i++) {
		vector<double> tArr(dataSize);
		for (int j = 0; j < dataSize; j++)
			tArr[j] = iArr[i][j];

		cout << "\nQuery:\t";
		for (int j = 0; j < tArr.size(); j++)
			cout << tArr[j] << "\t";
		cout << endl;

		auto rArr = get_empty_res(1, size);
		//rek_serial(gArr, size, depth, tArr, dataRows, dataSize, rArr, distances, resID);
		rek_amp(gArr, size, depth, tArr, dataRows, dataSize, rArr, distances, resID);
		double q = rek_quali(rArr, 0, size);
		sumQuali += q;
		double res = get_max_value_index(rArr, 0, size)*step;
		cout << "Result " << i << ": " << res << "/" << iArr[i][elements] << " Quali: " << q << endl;
		if (res != iArr[i][elements])
			errors++;
		//visualize(&rArr[0], size);
	}

	sumQuali /= dataRows;
	clock_t end = clock();

	cout << "GymyPoints:\t" << size << ", " << depth << endl;
	cout << "classes:\t" << classes << endl;
	cout << "samples:\t" << classes*samplesPerClass << endl;
	cout << "elements:\t" << elements << endl;
	cout << "elementsToIgnore:\t" << elem_to_ignore << endl;
	cout << "Average Rek Quali:\t" << sumQuali << endl;
	cout << "Rek Errors:\t" << errors << endl;
	cout << "Time:\t" << double(end - begin) / CLOCKS_PER_SEC << endl;
	cout << "Fertig" << endl;

}