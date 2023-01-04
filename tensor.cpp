#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */
using namespace std;


void Tensor::init_random(float mean, float std){

    if(data){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }

    }else{
        throw(tensor_not_initialized());
    }
}


Tensor::Tensor(){
    data = nullptr;
    r = 0;
    c = 0;
    d = 0;
}


Tensor::Tensor(int r, int c, int d, float v){
    //assign parameters
    this->r = r;
    this->c = c;
    this->d = d;

    if (!data){ //creating tensor
        data = new float**[r];
        for (int i = 0; i<r; i++){
            data[i] = new float*[c];
            for (int j = 0; j<c; j++){
                data[i][j] = new float[d];
                for(int k = 0; k<d; k++){
                    this->data[i][j][k] = v;
                }
            }
        }
    }else{

        //deleting older tensor
        for (int i = 0; i<r; i++){
            for (int j = 0; j<c; j++){
                delete[]data[i][j];
            }
        }

        for (int i = 0; i<r; i++){
            delete[] data[i];
        }

        delete[] data;

        //creating new tensor
        data = new float**[r];
        for (int i = 0; i<r; i++){
            data[i] = new float*[c];
            for (int j = 0; j<c; j++){
                data[i][j] = new float[d];
                for(int k = 0; k<d; k++){
                    this->data[i][j][k] = v;
                }
            }
        }

    }

}


Tensor::~Tensor(){

    if (data){

        //deleting tensor

        for (int i = 0; i<r; i++){
            for (int j = 0; j<c; j++){
                delete[]data[i][j];
            }
        }

        for (int i = 0; i<r; i++){
            delete[] data[i];
        }


        delete[] data;

        data = nullptr;
    }

}


float Tensor::operator()(int i, int j, int k) const{

    if(i > r || j > c || k > d){
        throw index_out_of_bound();
    }

    return data[i][j][k];
}


float& Tensor::operator()(int i, int j, int k){

    if(i > r || j > c || k > d){
        throw index_out_of_bound();
    }

    return data[i][j][k];
}


Tensor::Tensor(const Tensor& that){

    //assign values
    r = that.r;
    c = that.c;
    d = that.d;

    //creating Tensor
    data = new float**[that.r];
    for (int i = 0; i<that.r; i++){
        data[i] = new float*[that.c];
        for (int j = 0; j<that.c; j++){
            data[i][j] = new float[that.d];
            for (int k= 0; k<that.d; k++){
                data[i][j][k] = that.data[i][j][k];
            }
        }
    }

}


bool Tensor::operator==(const Tensor& rhs) const{

    if(r != rhs.r || c != rhs.c || d != rhs.d){
        throw dimension_mismatch();
    }else{

        bool result = true;
            for(int i = 0; i<r; i++ ){
                for (int j = 0; j<c; j++){
                    for (int k = 0; k<d; k++){
                        result = abs(data[i][j][k] - rhs.data[i][j][k]) < EPSILON ? result : false;
                    }
                }
            }

        return result;
    }
}


Tensor Tensor::operator-(const Tensor &rhs){

    if(r != rhs.r || c != rhs.c || d != rhs.d){
        throw dimension_mismatch();
    }else{

        Tensor result(this->r, this->c, this->d);
            for(int i = 0; i<rhs.r; i++){
                for(int j = 0; j<rhs.c; j++){
                    for(int k = 0; k<rhs.d; k++){
                        result.data[i][j][k] = data[i][j][k] - rhs.data[i][j][k];
                    }
                }
            }

        return result;
    }
}


Tensor Tensor::operator+(const Tensor &rhs){

    if(r != rhs.r || c != rhs.c || d != rhs.d){
            throw dimension_mismatch();
     }else{

        Tensor result(this->r, this->c, this->d, 0);
            for(int i = 0; i<rhs.r; i++){
                for(int j = 0; j<rhs.c; j++){
                    for(int k = 0; k<rhs.d; k++){
                        result.data[i][j][k] = data[i][j][k] + rhs.data[i][j][k];
                    }
                }
            }

        return result;
    }
}


Tensor Tensor::operator *(const Tensor &rhs){

    if(r != rhs.r || c != rhs.c || d != rhs.d){
            throw dimension_mismatch();
    }else{

        Tensor result(this->r, this->c, this->d);
            for(int i = 0; i<rhs.r; i++){
                for(int j = 0; j<rhs.c; j++){
                    for(int k = 0; k<rhs.d; k++){
                     result.data[i][j][k] = data[i][j][k] * rhs.data[i][j][k];
                    }
                }
            }

        return result;
    }
}


Tensor Tensor::operator/(const Tensor &rhs){

    if(r != rhs.r || c != rhs.c || d != rhs.d){
       throw dimension_mismatch();
    }else{

        Tensor result(this->r, this->c, this->d);
            for(int i = 0; i<rhs.r; i++){
                for(int j = 0; j<rhs.c; j++){
                    for(int k = 0; k<rhs.d; k++){
                        result.data[i][j][k] = data[i][j][k] / rhs.data[i][j][k];
                    }
                }
            }

        return result;
    }
}


Tensor Tensor::operator-(const float &rhs){

    Tensor result(this->r, this->c, this->d);
        for(int i = 0; i<this->r; i++){
            for(int j = 0; j<this->c; j++){
                for(int k = 0; k<this->d; k++){
                    result.data[i][j][k] = data[i][j][k] - rhs;
                }
            }
        }

    return result;
}


Tensor Tensor::operator+(const float &rhs){

    Tensor result(this->r, this->c, this->d);
        for(int i = 0; i<this->r; i++){
            for(int j = 0; j<this->c; j++){
                for(int k = 0; k<this->d; k++){
                    result.data[i][j][k] = data[i][j][k] + rhs;
                }
            }
        }

    return result;
}


Tensor Tensor::operator*(const float &rhs){

    Tensor result(r, c, d);
        for(int i = 0; i<r; i++){
            for(int j = 0; j<c; j++){
                for(int k = 0; k<d; k++){
                    result.data[i][j][k] = data[i][j][k] * rhs;
                }
            }
        }

    return result;
}


Tensor Tensor::operator/(const float &rhs){

    Tensor result(this->r, this->c, this->d);
        for(int i = 0; i<this->r; i++){
            for(int j = 0; j<this->c; j++){
                for(int k = 0; k<this->d; k++){
                    result.data[i][j][k] = data[i][j][k] / rhs;
                }
            }
        }

    return result;
}


Tensor& Tensor::operator=(const Tensor &other){

    if (data){

        //deleting older Tensor
        for (int i = 0; i<r; i++){
            for (int j = 0; j<c; j++){
                delete[]data[i][j];
            }
        }

        for (int i = 0; i<r; i++){
            delete[] data[i];
        }

        delete[] data;

    }

    //assign new values
     r = other.rows();
     c = other.cols();
     d = other.depth();

     //assign new Tensor
     data = new float**[r];
         for (int i = 0; i<r; i++){
             data[i] = new float*[c];
             for (int j = 0; j<c; j++){
                 data[i][j] = new float[d];
                 for(int k = 0; k<d; k++){
                     data[i][j][k] = other.data[i][j][k];
                 }
             }
         }

        return *this;
}


void Tensor::init(int r, int c, int d, float v){

    //creating new Tensor
    data = new float**[r];
    for (int i = 0; i<r; i++){
        data[i] = new float*[c];
        for (int j = 0; j<c; j++){
            data[i][j] = new float[d];
            for(int k = 0; k<d; k++){
                this->data[i][j][k] = v; //assign value v
            }
        }
    }
}


void Tensor::clamp(float low, float high){


    float min = this -> data[0][0][0];
    float max = this -> data[0][0][0];


    for(int i = 0; i<r; i++){ //computing min and max Tensor values
        for(int j = 0; j<c; j++){
            for(int k = 0; k<d; k++){
                min = data[i][j][k] < min ? data[i][j][k] : min;
                max = data[i][j][k] > max ? data[i][j][k] : max;
            }
        }
    }


    //clamp to low
    if (min < low){

        for(int i = 0; i<r; i++){
            for(int j = 0; j<c; j++){
                for(int k = 0; k<d; k++){
                    if(data[i][j][k] < low){
                        data[i][j][k] = low;
                    }
                }
            }
        }

    }else{

        for(int i = 0; i<r; i++){
            for(int j = 0; j<c; j++){
                for(int k = 0; k<d; k++){
                    if(data[i][j][k] == min){
                        data[i][j][k] = low;
                    }
                }
            }
        }

    }


    //clamp to high
    if (max > high){

        for(int i = 0; i<r; i++){
            for(int j = 0; j<c; j++){
                for(int k = 0; k<d; k++){
                    if(data[i][j][k] > high){
                        data[i][j][k] = high;
                    }
                }
            }
        }

    }else{

        for(int i = 0; i<r; i++){
            for(int j = 0; j<c; j++){
                for(int k = 0; k<d; k++){
                    if(data[i][j][k] == max){
                        data[i][j][k] = high;
                    }
                }
            }
        }

    }

}


void Tensor::rescale(float new_max){

    float min_k = data[0][0][0];
    float max_k = data[0][0][0];

    for(int i = 0; i<r; i++){
        for(int j = 0; j<c; j++){
            for(int k = 0; k<d; k++){
                min_k = min_k < data[i][j][k] ? min_k : data[i][j][k];
                max_k = max_k > data[i][j][k] ? max_k : data[i][j][k];
            }
        }
    }

    for(int i = 0; i<r; i++){
        for(int j = 0; j<c; j++){
            for(int k = 0; k<d; k++){
                data[i][j][k] = ((data[i][j][k] - min_k)/(max_k - min_k))*new_max;
            }
        }
    }

}


Tensor Tensor::padding(int pad_h, int pad_w){

   Tensor new_Tensor = Tensor(r+(2*pad_h), c+(2*pad_w), d, 0.0);

    for(int k = 0; k<d; k++){
        for(int i = 0; i<r; i++){
            for(int j = 0; j<c; j++){
                new_Tensor.data[i+pad_h][j+pad_w][k] = data[i][j][k];

            }
        }
    }

    return new_Tensor;
}


Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end){

    Tensor new_Tensor(row_end - row_start, col_end - col_start, depth_end - depth_start, 0);

    for(int i = row_start; i<row_end; i++){
        for(int j = col_start; j<col_end; j++){
            for(int k = depth_start; k<depth_end; k++){
                new_Tensor.data[i - row_start][j - col_start][k - depth_start] = data[i][j][k];
            }
        }
    }

    return new_Tensor;
}


Tensor Tensor::concat(const Tensor &rhs, int axis){

    switch(axis){

        case 0:

            if ((c != rhs.c) || (d != rhs.d)){

                throw concat_wrong_dimension();

            }else{

                Tensor result( r+rhs.r , c, d, 0.0);

                for(int i = 0; i<r; i++){
                    for(int j = 0; j<c; j++){
                        for(int k = 0; k<d; k++){
                            result.data[i][j][k] = data[i][j][k];
                        }
                    }
                }

                for(int i = 0; i<rhs.r; i++){
                    for(int j = 0; j<rhs.c; j++){
                        for(int k = 0; k<rhs.d; k++){
                            result.data[r + i][j][k] = rhs.data[i][j][k];
                        }
                    }
                }

                return result;
            }
        case 1:

            if ((r != rhs.r) || (d != rhs.d)){

                throw concat_wrong_dimension();

            }else{

                Tensor result(r, c+rhs.c, d, 0.0);

                for(int i = 0; i<r; i++){
                    for(int j = 0; j<c; j++){
                        for(int k = 0; k<d; k++){
                            result.data[i][j][k] = data[i][j][k];
                        }
                    }
                }

                for(int i = 0; i<rhs.r; i++){
                    for(int j = 0; j<rhs.c; j++){
                        for(int k = 0; k<rhs.d; k++){
                            result.data[i][c + j][k] = rhs.data[i][j][k];
                        }
                    }
                }

                return result;
            }

        case 2:

            if ((r != rhs.r) || (d != rhs.d)){

                throw concat_wrong_dimension();

            }else{

                Tensor result( r, c, d + rhs.d, 0.0);

                for(int i = 0; i<r; i++){
                    for(int j = 0; j<c; j++){
                        for(int k = 0; k<d; k++){
                            result.data[i][j][k] = data[i][j][k];
                        }
                    }
                }

                for(int i = 0; i<rhs.r; i++){
                    for(int j = 0; j<rhs.c; j++){
                        for(int k = 0; k<rhs.d; k++){
                            result.data[i][j][d + k] = rhs.data[i][j][k];
                        }
                    }
                }

                return result;
            }

        default:
            throw concat_wrong_dimension();

    }

}


float Tensor::_convolve_aux(const Tensor &part, const Tensor &filter, int r, int c, int layer){

    float result = 0;

        for(int i = r; i<filter.r + r; i++){
            for (int j = c; j<filter.c + c; j++){
                result += filter.data[i - r][j - c][layer] * part.data[i][j][layer];
            }
        }

        return result;
}


Tensor Tensor::convolve(const Tensor &f){

    if((f.c*f.r)%2!=0 && d==f.d){

        Tensor result(r, c, d, 0.0);
        Tensor padded;
        padded = this->padding((f.r - 1)/2, (f.c - 1)/2);

        for(int k = 0; k<d; k++){
            for(int i = 0; i<r; i++){
                for(int j = 0; j<c; j++){
                    result.data[i][j][k] = _convolve_aux(padded, f, i, j, k);
                }
            }
        }

        return result;

    }else{

        throw filter_odd_dimensions();

    }

}


/* UTILITY */


int Tensor::rows() const{
    return r;
}


int Tensor::cols()const{
    return c;
}


int Tensor::depth()const{
    return d;
}


float Tensor::getMin(int k)const{

    float min_k = data[0][0][k];

        for( int i = 0; i<r; i++ ){
            for (int j = 0; j<c; j++){
                min_k = data[i][j][k] < min_k ? data[i][j][k] : min_k;
            }
        }

    return min_k;
}


float Tensor::getMax(int k)const{

    float max_k = data[0][0][k];

        for( int i = 0; i<r; i++ ){
            for (int j = 0; j<c; j++){
                max_k = data[i][j][k] > max_k ? data[i][j][k] : max_k;
            }
        }

    return max_k;

}


void Tensor::showSize() const{
    cout<<r<<" x "<<c<<" x "<<d;
}

/* IOSTREAM */


ostream& operator<<(ostream& stream, const Tensor& obj){

        for(int i = 0; i<obj.r; i++ ){

            for (int k = 0; k<obj.d; k++){

                stream<<"[ ";

                    for (int j = 0; j<obj.c; j++){
                        stream<<obj.data[i][j][k]<<" ";
                    }

                stream<<"]";
                stream<<"  ";

            }

            stream<<endl;

        }

    return stream;

}


void Tensor::read_file(string filename){

    ifstream file;
    file.open(filename);


    if(!file){

        throw unable_to_read_file();

    }else{

        while(!file.eof()){


            file>>r;
            file>>c;
            file>>d;


            if (data){ //deleting older tensor

                for (int i = 0; i<r; i++){
                    for (int j = 0; j<c; j++){
                        delete[]data[i][j];
                    }
                }

                for (int i = 0; i<r; i++){
                    delete[] data[i];
                }

                delete[] data;

                data = nullptr;

            }


            if (!data){ //assign new tensor

                data = new float**[r];
                for (int i = 0; i<r; i++){
                    data[i] = new float*[c];
                    for (int j = 0; j<c; j++){
                        data[i][j] = new float[d];
                    }
                }

            }


            //assign values to the Tensor
            for(int i = 0; i<r; i++){
                for (int j = 0; j<c; j++){
                    for (int k = 0; k<d; k++){
                        file>>data[i][j][k];
                    }
                }
            }


        }

    }

    file.close();
}


void Tensor::write_file(string filename){

    ofstream file;
    file.open(filename);


    if(!file){

        throw unable_to_read_file();

    }else {


        file << r << endl;
        file << c << endl;
        file << d << endl;

        for(int i = 0; i<r; i++ ){
            for (int j = 0; j<c; j++){
                for (int k = 0; k<d; k++){
                    file<<data[i][j][k]<<endl;
                }
            }
        }


    }

    file.close();
}
