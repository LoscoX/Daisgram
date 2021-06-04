#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

/* TEMPLATE PER SCORRERE L'ARRAY

for(int i = 0; i<r; i++ ){
    for (int j = 0; j<c; j++){
        for (int k = 0; k<d; k++){
            data[i][j][k];
        }
    }
} */




/**
 * Load a bitmap from file
 *
 * @param filename String containing the path of the file
 */

void DAISGram::load_image(string filename){
    BmpImg img = BmpImg();

    img.read(filename.c_str());

    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for(int i=0;i<img.get_height();i++){
        for(int j=0;j<img.get_width();j++){
            data(i,j,0) = (float) img.red_at(j,i);
            data(i,j,1) = (float) img.green_at(j,i);
            data(i,j,2) = (float) img.blue_at(j,i);
        }
    }
}


/**
 * Save a DAISGram object to a bitmap file.
 *
 * Data is clamped to 0,255 before saving it.
 *
 * @param filename String containing the path where to store the image.
 */
void DAISGram::save_image(string filename){

    data.clamp(0,255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for(int i=0;i<getRows();i++){
        for(int j=0;j<getCols();j++){
            img.set_pixel(j,i,(unsigned char) data(i,j,0),(unsigned char) data(i,j,1),(unsigned char) data(i,j,2));
        }
    }

    img.write(filename);

}


/**
 * Generate Random Image
 *
 * Generate a random image from nois
 *
 * @param h height of the image
 * @param w width of the image
 * @param d number of channels
 * @return returns a new DAISGram containing the generated image.
 */
void DAISGram::generate_random(int h, int w, int d){
    data = Tensor(h,w,d,0.0);
    data.init_random(128,50);
    data.rescale(255);
}

        DAISGram::DAISGram(){
        }

        DAISGram::~DAISGram(){
        }

        /**
         * Get rows
         *
         * @return returns the number of rows in the image
         */
        int DAISGram::getRows(){
            return data.rows();
        }

        /**
         * Get columns
         *
         * @return returns the number of columns in the image
         */
        int DAISGram::getCols(){
            return data.cols();
        }

        /**
         * Get depth
         *
         * @return returns the number of channels in the image
         */
        int DAISGram::getDepth(){
            return data.depth();
        }

        /**
         * Brighten the image
         *
         * It sums the bright variable to all the values in the image.
         *
         * Before returning the image, the corresponding tensor should be clamped in [0,255]
         *
         * @param bright the amount of bright to add (if negative the image gets darker)
         * @return returns a new DAISGram containing the modified object
         */
        DAISGram DAISGram::brighten(float bright){
            DAISGram new_tensor;
            new_tensor.data = data;
            for(int i = 0; i<data.rows(); i++ ){
                for (int j = 0; j<data.cols(); j++){
                    for (int k = 0; k<data.depth(); k++){
                         new_tensor.data(i,j,k)= data(i,j,k) + bright;
                    }
                }
            }
            return new_tensor;
        }

        /**
         * Create a grayscale version of the object
         *
         * A grayscale image is produced by substituting each pixel with its average on all the channel
         *
         * @return returns a new DAISGram containing the modified object
         */
        DAISGram DAISGram::grayscale(){
            DAISGram result;

           result.data = data;

            for (int i = 0; i <data.rows(); i++){
                for (int j = 0; j<data.cols(); j++){
                    float average;
                    for (int k = 0; k<data.depth(); k++){
                        average += result.data(i,j,k);
                    }
                    average /= data.depth();
                    for (int k = 0; k<data.depth(); k++){
                        result.data(i,j,k) = average;
                    }
                }
            }

            return result;
        }

        void DAISGram::swap(int k1, int k2){
            DAISGram supp;
            supp.data = Tensor(data.rows(), data.cols(), 0, 0.0);
            for(int i = 0; i<data.rows(); i++ ){
                for (int j = 0; j<data.cols(); j++){
                    supp.data(i,j,0)= data(i,j,k1);
                }
            }
             for(int i = 0; i<data.rows(); i++ ){
                for (int j = 0; j<data.cols(); j++){
                    data(i,j,k1) = data(i,j,k2);
                }
            }
             for(int i = 0; i<data.rows(); i++ ){
                for (int j = 0; j<data.cols(); j++){
                    data(i,j,k2) = supp.data(i,j,0);
                }
            }
        }


        /**
         * Create a Warhol effect on the image
         *
         * This function returns a composition of 4 different images in which the:
         * - top left is the original image
         * - top right is the original image in which the Red and Green channel are swapped
         * - bottom left is the original image in which the Blue and Green channel are swapped
         * - bottom right is the original image in which the Red and Blue channel are swapped
         *
         * The output image is twice the dimensions of the original one.
         *
         * @return returns a new DAISGram containing the modified object
         */
        DAISGram DAISGram::warhol(){
            DAISGram new_tensor;
            new_tensor.data = Tensor(data.rows()*2, data.cols()*2, data.depth()*2, 0.0);
            for(int i = 0; i<data.rows(); i++){                                       //in alto a sinistra
                for (int j = 0; j<data.cols(); j++){
                    for (int k = 0; k<data.depth(); k++){
                        new_tensor.data(i,j,k) = data(i,j,k);
                    }
                }
            }
            swap(0,1);
            for(int i = 0; i<data.rows(); i++){                           //in alto a destra
                for (int j = data.cols() ; j<data.cols()*2; j++){
                    for (int k = 0; k<data.depth(); k++){
                        new_tensor.data(i,j,k) = data(i, j - data.cols(),k);
                    }
                }
            }
            swap(0,1);
            swap(1,2);
            for(int i = data.rows(); i<data.rows()*2; i++){                                      //in basso a sinistra
                for (int j = 0; j<data.cols(); j++){
                    for (int k = 0; k<data.depth(); k++){
                        new_tensor.data(i,j,k) = data( i - data.rows(), j, k);
                    }
                }
            }
            swap(1,2);
            swap(0,2);
            for(int i = data.rows(); i<data.rows()*2; i++){                          //in basso a destra
                for (int j = data.cols(); j<data.cols()*2; j++){
                    for (int k = 0; k<data.depth(); k++){
                        new_tensor.data(i,j,k) = data(i - data.rows(), j - data.cols(),k);
                    }
                }
            }

            return new_tensor;
        }

        /**
         * Sharpen the image
         *
         * This function makes the image sharper by convolving it with a sharp filter
         *
         * filter[3][3]
         *    0  -1  0
         *    -1  5 -1
         *    0  -1  0
         *
         * Before returning the image, the corresponding tensor should be clamped in [0,255]
         *
         * @return returns a new DAISGram containing the modified object
         */
         /*
        DAISGram DAISGram::sharpen(){
            DAISGram new_tensor;
            new_tensor.data = Tensor(data.rows(), data.cols(), data.depth(), 0.0);
			DAISGram filter;
            filter.data = Tensor(3,3,1, 0.0);
            for(int k=0; k<data.depth(); k++){
              filter.data(0,0,k) = 0.0;
              filter.data(0,1,k) = -1.0;
              filter.data(1,0,k) = -1.0;
              filter.data(1,1,k) = 5.0;
              filter.data(2,1,k) = -1.0;
              filter.data(1,2,k) = -1.0;
              filter.data(0,2,k) = 0.0;
              filter.data(2,0,k) = 0.0;
              filter.data(2,2,k) = 0.0;
            }
            new_tensor.data = data.convolve(filter.data);
			new_tensor.data.clamp(0,255);
            return new_tensor;
        }
        */
        DAISGram DAISGram::sharpen(){
            DAISGram new_tensor;
            new_tensor.data = Tensor(data.rows(), data.cols(), data.depth(), 0.0);
			DAISGram filter;
             filter.data = Tensor(3,3,data.depth(), 0.0);
             for(int k=0; k<data.depth(); k++){
                 filter.data(0,0,k) = 0;
                 filter.data(0,1,k) = -1.0;
                 filter.data(1,0,k) = -1.0;
                 filter.data(1,1,k) = 5.0;
                 filter.data(2,1,k) = -1.0;
                 filter.data(1,2,k) = -1.0;
                 filter.data(0,2,k) = 0.0;
                 filter.data(2,0,k) = 0.0;
                 filter.data(2,2,k) = 0.0;
             }
            new_tensor.data = data.convolve(filter.data);
            new_tensor.data.clamp(0,255);
            new_tensor.data.rescale(255);
            return new_tensor;
        }

        /**
         * Emboss the image
         *
         * This function makes the image embossed (a light 3D effect) by convolving it with an
         * embossing filter
         *
         * filter[3][3]
         *    -2 -1  0
         *    -1  1  1
         *     0  1  2
         *
         * Before returning the image, the corresponding tensor should be clamped in [0,255]
         *
         * @return returns a new DAISGram containing the modified object
         */
        DAISGram DAISGram::emboss(){
            DAISGram new_tensor;
            new_tensor.data = Tensor(data.rows(), data.cols(), data.depth(), 0.0);
            DAISGram filter;
            filter.data = Tensor(3,3,data.depth(), 0.0);
            for(int k=0; k<data.depth(); k++){
              filter.data(0,0,k) = -2.0;
              filter.data(0,1,k) = -1.0;
              filter.data(1,0,k) = -1.0;
              filter.data(1,1,k) = 1.0;
              filter.data(2,1,k) = 1.0;
              filter.data(1,2,k) = 1.0;
              filter.data(0,2,k) = 0.0;
              filter.data(2,0,k) = 0.0;
              filter.data(2,2,k) = 2.0;
            }
            new_tensor.data = data.convolve(filter.data);
            new_tensor.data.clamp(0,255);
            new_tensor.data.rescale(255);
            return new_tensor;
        }

        /**
         * Smooth the image
         *
         * This function remove the noise in an image using convolution and an average filter
         * of size h*h:
         *
         * c = 1/(h*h)
         *
         * filter[3][3]
         *    c c c
         *    c c c
         *    c c c
         *
         * @param h the size of the filter
         * @return returns a new DAISGram containing the  object
         */
        DAISGram DAISGram::smooth(int h){
            DAISGram new_tensor;
            new_tensor.data = data;
            DAISGram filter;
            filter.data = Tensor(3,3,data.depth(), 0.0);
            float c = 1.0 / (float) (h*h);
            for(int k=0; k<data.depth(); k++){
              filter.data(0,0,k) = c;
              filter.data(0,1,k) = c;
              filter.data(1,0,k) = c;
              filter.data(1,1,k) = c;
              filter.data(2,1,k) = c;
              filter.data(1,2,k) = c;
              filter.data(0,2,k) = c;
              filter.data(2,0,k) = c;
              filter.data(2,2,k) = c;
            }
            new_tensor.data = data.convolve(filter.data);
            new_tensor.data.rescale(255);
            new_tensor.data.clamp(0,255);
            return new_tensor;
        }

        /**
         * Edges of an image
         *
         * This function extract the edges of an image by using the convolution
         * operator and the following filter
         *
         *
         * filter[3][3]
         * -1  -1  -1
         * -1   8  -1
         * -1  -1  -1
         *
         * Remeber to convert the image to grayscale before running the convolution.
         *
         * Before returning the image, the corresponding tensor should be clamped in [0,255]
         *
         * @return returns a new DAISGram containing the modified object
         */
        DAISGram DAISGram::edge(){
            DAISGram new_tensor;
            new_tensor.data = Tensor(data.rows(), data.cols(), data.depth(), 0.0);
            DAISGram filter;
            filter.data = Tensor(3,3,data.depth(), 0.0);
            for(int k=0; k<data.depth(); k++){
              filter.data(0,0,k) = -1.0;
              filter.data(0,1,k) = -1.0;
              filter.data(1,0,k) = -1.0;
              filter.data(1,1,k) = 8.0;
              filter.data(2,1,k) = -1.0;
              filter.data(1,2,k) = -1.0;
              filter.data(0,2,k) = -1.0;
              filter.data(2,0,k) = -1.0;
              filter.data(2,2,k) = -1.0;
            }
            data = grayscale().data;
            new_tensor.data = data.convolve(filter.data);
            new_tensor.data.clamp(0,255);
            new_tensor.data.rescale(255);
            return new_tensor;
        }

        /**
         * Blend with anoter image
         *
         * This function generate a new DAISGram which is the composition
         * of the object and another DAISGram object
         *
         * The composition follows this convex combination:
         * results = alpha*this + (1-alpha)*rhs
         *
         * rhs and this obejct MUST have the same dimensions.
         *
         * @param rhs The second image involved in the blending
         * @param alpha The parameter of the convex combination
         * @return returns a new DAISGram containing the blending of the two images.
         */
        DAISGram DAISGram::blend(const DAISGram & rhs, float alpha){
            DAISGram new_tensor;
            new_tensor.data = Tensor(rhs.data.rows(), rhs.data.cols(), rhs.data.depth(), 0.0);
                for(int i = 0; i<rhs.data.rows(); i++ ){
                    for (int j = 0; j<rhs.data.cols(); j++){
                        for (int k = 0; k<rhs.data.depth(); k++){
                            new_tensor.data(i,j,k) = rhs.data(i,j,k);
                        }
                    }
                }



            if(data.rows() == rhs.data.rows() && data.cols() == rhs.data.cols() && data.depth() == rhs.data.depth()){
              for(int i = 0; i<rhs.data.rows(); i++ ){
                  for (int j = 0; j<rhs.data.cols(); j++){
                      for (int k = 0; k<rhs.data.depth(); k++){
                          new_tensor.data(i,j,k) = ((data(i,j,k)*alpha)+(new_tensor.data(i,j,k))*(1-alpha));
                      }
                  }
              }

            }
            return new_tensor;
        }

        /**
         * Green Screen
         *
         * This function substitutes a pixel with the corresponding one in a background image
         * if its colors are in the surrounding (+- threshold) of a given color (rgb).
         *
         * (rgb - threshold) <= pixel <= (rgb + threshold)
         *
         *
         * @param bkg The second image used as background
         * @param rgb[] The color to substitute (rgb[0] = RED, rgb[1]=GREEN, rgb[2]=BLUE)
         * @param threshold[] The threshold to add/remove for each color (threshold[0] = RED, threshold[1]=GREEN, threshold[2]=BLUE)
         * @return returns a new DAISGram containing the result.
         */
        DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[]){
            DAISGram result;
            result.data = data;

            for (int i = 0; i<result.data.rows(); i++){
                for (int j = 0; j<result.data.cols(); j++){
                    bool in_threshold = true;
                    for (int k = 0; k<result.data.depth(); k++){
                        in_threshold = result.data(i, j, k) <= rgb[k]+threshold[k] ? in_threshold : false;
                        in_threshold = result.data(i, j, k) >= rgb[k]-threshold[k] ? in_threshold : false;
                    }
                    if (in_threshold){
                        for (int k = 0; k < result.data.depth(); k++){
                            result.data(i, j, k) = bkg.data(i, j, k);
                        }
                    }
                }
            }
            return result;
        }

        /**
         * Equalize
         *
         * Stretch the distribution of colors of the image in order to use the full range of intesities.
         *
         * See https://it.wikipedia.org/wiki/Equalizzazione_dell%27istogramma
         *
         * @return returns a new DAISGram containing the equalized image.
         */
        DAISGram DAISGram::equalize(){
            /*DAISGram result;
            result.data = data;
            for (int k= 0; k<data.depth(); k++){
                int intensity[256];
                for (int i = 0; i<256; i++){
                    intensity[i] = 0;
                }

                for (int i = 0; i<data.rows(); i++){
                    for (int j = 0; j<data.cols(); j++){
                        ++intensity[data(i, j, k)];
                    }
                }
                
                
            }*/

        }
