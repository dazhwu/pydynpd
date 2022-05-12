//
// Created by Tiger on 5/4/2022.
//

#include "Common_Functions.h"
#define EIGEN_INITIALIZE_MATRICES_BY_NAN

using namespace std;
using namespace Eigen;

vector<bool> row_has_nan(const RowMatrixXd &x)
{
    int num_rows=x.rows();
    vector<bool> tbr(num_rows);
    for (int i=0; i<num_rows; ++i){

        tbr[i]=!((x.row(i).array()==x.row(i).array()).all());
        
    }
   return tbr;
}

RowMatrixXd common_inv(Ref <RowMatrixXd> ori) {

    Eigen::FullPivLU <RowMatrixXd> Lu(ori);
	return ori.completeOrthogonalDecomposition().pseudoInverse();
//    if (Lu.isInvertible())
//        return ori.inverse();
//    else
//        return ori.completeOrthogonalDecomposition().pseudoInverse();


}

double standard_normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x / std::sqrt(2)) / 2;
}

/* The function checks if the array elements are consecutive
If elements are consecutive, then returns true, else returns
false */
bool areConsecutive(int* arr, int n)  //int arr[]
{
	//int n = sizeof(arr) / sizeof(arr[0]);
	if (n <  1)
		return false;

	int min = *min_element(arr, arr + n);

	int max = *max_element(arr, arr + n);

	if (max - min + 1 == n)
	{
		/* Create a temp array to hold visited flag of all elements.
		Note that, calloc is used here so that all values are initialized
		as false */
		bool *visited = (bool *)calloc(n, sizeof(bool));
		int i;
		for (i = 0; i < n; i++)
		{
			/* If we see an element again, then return false */
			if (visited[arr[i] - min] != false)
				return false;

			/* If visited first time, then mark the element as visited */
			visited[arr[i] - min] = true;
		}

		/* If all elements occur once, then return true */
		return true;
	}

	return false; // if (max - min  + 1 != n)
}

 
//void tokenize(std::string const &str, const char* delim,
//            std::vector<std::string> &out)
//{
//    char *token = strtok_s(const_cast<char*>(str.c_str()), delim);
//    while (token != nullptr)
//    {
//        out.push_back(std::string(token));
//        token = strtok_s(nullptr, delim);
//    }
//}

vector<string> splitString(string str, char splitter) {
	vector<string> result;
	string current = "";
	for (int i = 0; i < str.size(); i++) {
		if (str[i] == splitter) {
			if (current != "") {
				result.push_back(current);
				current = "";
			}
			continue;
		}
		current += str[i];
	}
	if (current.size() != 0)
		result.push_back(current);
	return result;
}
bool variable_exists(string var_name,  vector<string> &cols) {
	//for (string var_name : names)
	//	if (var_name == var_name)
	//		return true;

	//return false;


	// if (std::find(cols.begin(), cols.end(), var_name) == cols.end())
	// 	return false;
	// return true;
	if (getIndex(cols, var_name)>=0)
		return true;
	else
		return false;

}

int getIndex(vector<string> &v, string K)
{
    auto it = std::find(v.begin(), v.end(), K);
 
    // If element was found
    if (it != v.end())
    {
     

        int index = it - v.begin();
        return index;
    }
    else {
        // If the element is not
        // present in the vector
        return -1;
    }
}

void lag(Ref<RowMatrixXd> mat, Ref<RowMatrixXd> lagged, int N, int lag_number, double fill){
	int	height = mat.rows() / N;
	int width = mat.cols();
	//int start_row, end_row;

	for (int i=0;i<N; ++i){
		//start_row = i * height;
		//end_row = start_row + height;
		Ref<RowMatrixXd> mat_i = mat.block(i*height, 0, height, width);
		Ref<RowMatrixXd> lagged_i = lagged.block(i*height, 0, height, width);

		if (!isnan(fill))
			lagged_i.block(0,0, lag_number, width) = RowMatrixXd::Zero(lag_number, width);
		else
			lagged_i.block(0,0, lag_number, width) = RowMatrixXd(lag_number, width);

		lagged_i.block(lag_number,0, height-lag_number, width) = mat_i.block(0,0,height - lag_number, width);
	}
}

RowMatrixXd get_first_diff_table(Ref<RowMatrixXd> ori_arr, int N){
	int num_cols = ori_arr.cols();
	int num_rows = ori_arr.rows();
	//int height = num_rows / N;

	RowMatrixXd lag_arr (num_rows, num_cols);
	RowMatrixXd tbr_arr (num_rows, num_cols);

	lag(ori_arr, lag_arr, N, 1, NAN);

	tbr_arr = ori_arr - lag_arr;
	return tbr_arr;

}



RowMatrixXd	get_fod_table(Ref<RowMatrixXd> ori_arr, int N){
	int num_cols = ori_arr.cols();
	int num_rows = ori_arr.rows();
	int height = num_rows / N;

	RowMatrixXd tbr (num_rows, num_cols);
	RowMatrixXd next_sum (1, num_cols);
	RowMatrixXd this_sum (1, num_cols);
	RowMatrixXd this_avg (1, num_cols);
	

//	tbr = np.empty((num_rows, num_cols), dtype='float64')
//	next_sum = np.empty((1, num_cols), dtype='float64')
//	this_sum = np.empty((1, num_cols), dtype='float64')
//	this_avg = np.empty((1, num_cols), dtype='float64')
//	temp = np.empty((height, num_cols), dtype='float64')

//	tbr[:] = np.NaN

//	this_sum[:] = np.NaN
	int this_count;
	for (int i=0; i<N; ++i){
		Ref<RowMatrixXd> ori_i = ori_arr.block(i * height,0, height, num_cols);
		Ref<RowMatrixXd> tbr_i= tbr.block(i * height,0, height, num_cols);
		
		RowMatrixXd temp (height, num_cols);
		RowMatrixXd next_sum (1, num_cols);
		int next_count = 0;
		//for j in range(height - 2, -1, -1):
		
		for (int j=height-2; j>=0; --j){

			if (!(ori_i.row(j).array()==ori_i.row(j).array()).all()){
				this_count = next_count;
				this_sum = next_sum;
				temp.row(j) = temp.row(j + 1);
			}else{
				this_count = next_count + 1;
				for(int k=0; k<num_cols; ++k){
					double temp_next_sum, temp_ori;
					if (isnan(next_sum(0,k)))
						temp_next_sum=0;
					else
						temp_next_sum=next_sum(0,k);
					
					if (isnan(ori_i(j+1,k)))
						temp_ori=0;
					else
						temp_ori=ori_i(j+1,k);

					this_sum(0,k)=temp_next_sum+temp_ori;
				}
				
				this_avg = this_sum * (1.0 / this_count);
				temp.row(j) = (ori_i.row(j) - this_avg) * sqrt(this_count / (this_count + 1));
			}
			next_sum = this_sum;
			next_count = this_count;

			
		}



	
	tbr_i.row(0) = RowMatrixXd(1, num_cols);
	tbr_i.block(1,0, height-1, num_cols) = temp.block(0, 0, height - 1, num_cols);

	}


	return tbr;
}

void saveData(string fileName, RowMatrixXd  matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
 
    ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}
