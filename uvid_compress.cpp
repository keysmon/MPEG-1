/* uvid_compress.cpp
   CSC 485B/578B - Data Compression - Summer 2023

   Starter code for Assignment 4

   Reads video data from stdin in uncompresed YCbCr (YUV) format 
   (With 4:2:0 subsampling). To produce this format from 
   arbitrary video data in a popular format, use the ffmpeg
   tool and a command like 

     ffmpeg -i videofile.mp4 -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./this_program <width> height>

   Note that since the width/height of each frame is not encoded into the raw
   video stream, those values must be provided to the program as arguments.

   B. Bird - 2023-07-08
*/
#include <cstdint>
#include <math.h> 
#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <tuple>
#include "output_stream.hpp"
#include "yuv_stream.hpp"
#include "uvg_common.hpp"
#include "huffman_tree.hpp"
#include "cmath"
#include <algorithm>

using namespace std;

vector<vector<unsigned char>> Y_buffer;
vector<vector<unsigned char>> Cb_buffer;
vector<vector<unsigned char>> Cr_buffer;
OutputBitStream output_stream {std::cout};
int counter = 0;
int total_frames = 0;
const int search_range = 3;
const int macro_size = 16;
const int block_size = 8;
vector<vector<int>> block_cords;
const double pi = 3.14159265358979323846;
int write_data = 1;
const int DCT_Tune = 3;

vector<vector<double>> quantization_matrix_lumin = {
{16,11,10,16,24,40,51,61},
{12,12,14,19,26,58,60,55},
{14,13,16,24,40,57,69,56},
{14,17,22,29,51,87,80,62},
{18,22,37,56,68,109,103,77},
{24,35,55,64,81,104,113,92},
{49,64,78,87,103,121,120,101},
{72,92,95,98,112,100,103,99}};

vector<vector<double>> quantization_matrix_chrom = {
{17,18,24,47,99,99,99,99},
{18,21,26,66,99,99,99,99},
{24,26,56,99,99,99,99,99},
{47,66,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99}};

vector<vector<int>> zig_zag_route = {
{0,0},
{0,1},{1,0},
{0,2},{1,1},{2,0},
{0,3},{2,1},{1,2},{3,0},
{0,4},{3,1},{2,2},{1,3},{4,0},
{0,5},{4,1},{3,2},{2,3},{1,4},{5,0},
{0,6},{5,1},{4,2},{3,3},{2,4},{1,5},{6,0},
{0,7},{6,1},{2,5},{3,4},{4,3},{5,2},{1,6},{7,0},
{1,7},{2,6},{3,5},{4,4},{5,3},{6,2},{7,1},
{2,7},{3,6},{4,5},{5,4},{6,3},{7,2},
{3,7},{4,6},{5,5},{6,4},{7,3},
{4,7},{5,6},{6,5},{7,4},
{5,7},{6,6},{7,5},
{6,7},{7,6},
{7,7}};

vector<vector<int>> zig_zag_coords(int height, int width,int frequent_size){
    

     vector<vector<int>> coords = {};
    int num_block_height = int(ceil(double(height)/8));
    int num_block_width = int(ceil(double(width)/8));
    
    for (int y = 0;y < num_block_height ; y++){
        for (int x = 0 ; x < num_block_width ; x++){
            for (int i = 0 ; i< frequent_size;i++){
                int y_cord = zig_zag_route[i][0] + y * block_size;
                int x_cord = zig_zag_route[i][1] + x * block_size;
                if (y_cord < height && x_cord < width){
                    vector<int> temp = {y_cord,x_cord};
                    coords.push_back(temp);
                }
            }
        }
    }
    return coords;

}

vector<vector<double>> dctTransform (vector<vector<double>> matrix){
    int x,y,u,v;
    int len = 8;
    // dct will store the discrete cosine transform
    vector<vector<double>> dct = create_2d_vector<double>(block_size,block_size);
    double au, av, dct1, sum;
    for (u = 0; u < len; u++) {
        for (v = 0; v < len; v++) {
            
            // ci and cj depends on frequency as well as
            // number of row and columns of specified matrix
            if (u == 0)
                au = 1/sqrt(2);
            else
                au = 1;
            if (v == 0)
                av = 1/sqrt(2);
            else
                av = 1;

            // sum will temporarily store the sum of
            // cosine signals
            sum = 0;
            for (x = 0; x < block_size; x++) {
                for (y = 0; y < block_size; y++) {
                    dct1 = matrix[x][y] *
                           double(cos((2 * x + 1) * u * pi / (2 * block_size))) *
                           double(cos((2 * y + 1) * v * pi / (2 * block_size)));
                    sum = sum + dct1;
                }
            }
            dct[u][v] = au * av * sum * 0.25;
        }
    }
    return dct;
}

vector<vector<double>> inverse_dctTransform(vector<vector<double>> matrix)
{   
    
    int u,v,x,y;
   
    // dct will store the discrete cosine transform
    vector<vector<double>> dct = create_2d_vector<double>(block_size,block_size);
    double au, av, dct1, sum;
    for (x = 0;  x < block_size; x++) {
        for (y = 0; y < block_size; y++) {
            
            sum = 0;
            for (u = 0; u < block_size; u++) {
                for (v = 0; v < block_size; v++) {
                    
                    if (u == 0)
                        au = 1/sqrt(2);
                    else
                        au = 1;
                    if (v == 0)
                        av = 1/sqrt(2);
                    else
                        av = 1;

                    dct1 = matrix[u][v] *
                           cos((2 * x + 1) * u * pi / (2 * block_size)) *
                           cos((2 * y + 1) * v * pi / (2 * block_size));
                    sum = sum + dct1 * au * av;
                }
            }
            dct[x][y] =  sum * 0.25;
        }
    }
    
    return dct;
}

vector<unsigned char> zig_zag_flatten (vector<vector<unsigned char>> image){
    vector<unsigned char>  res = {};
    for (auto pair:zig_zag_route){
        res.push_back(image[pair[0]][pair[1]]);
    }       
    return res;
} 

vector<unsigned char> delta_compress(vector<unsigned char> flattened){
    vector<unsigned char> result(flattened.size(), 0);
    int last = 0;
    result[0] = flattened[0];
    for (int i = 1; i < flattened.size(); i++)
    {
        unsigned char current = flattened[i];
        result[i] =  (flattened[i-1] - current)/2 + 128;
        
    }
    return result;
}

vector<unsigned char> delta_decompress(vector<unsigned char> input){
    vector<unsigned char> res = {};
    res.push_back(input[0]);
    for (int i=1;i<input.size();i++){
        res.push_back(res[i-1] - (input[i] - 128)*2);
    }

    return res;
}

vector<unsigned char> RLE_without_first_element (vector<vector<unsigned char>> image){
    int checked = 0;
    vector<unsigned char> flattened_in_zig_zag = {};
    vector<unsigned char> result = {};
    

    //for (auto pair:block_cords){
    for (int i=1;i<64;i++){  
        auto pair = zig_zag_route[i];
        int y_cord = pair[0];
        int x_cord = pair[1];
        flattened_in_zig_zag.push_back(image[y_cord][x_cord]);
    }



    for( int i = 0; i < flattened_in_zig_zag.size(); i++){
        int count = 1;
        //bool temp = (flattened_in_zig_zag[i] == 128 || flattened_in_zig_zag[i] == 0);
        while(flattened_in_zig_zag[i] == flattened_in_zig_zag[i + 1] && i < flattened_in_zig_zag.size() - 1 ){
            count++;
            i++;
        }
        
        result.push_back(flattened_in_zig_zag[i]);
        result.push_back(count);
    }
    return result;
} 

vector<unsigned char> huffman_compress (vector<unsigned char> image){
    vector<unsigned char> result = {};
    huf_tree test;
    image.push_back(255);
    
	string sss = "";
	int count = 0;
	unsigned char temp;
    
    
    for ( auto s:image){

        
        sss += test.addnode(s);
        while(sss.size() - count >= 8) {
        	for(int i = count; i < count + 8; i ++) {
        		temp <<= 1;
        		if(sss[i] == '1') temp |= 1;
			}
            result.push_back(temp);
			count += 8;
		}
	}

    int left_len = sss.size() - count;
    if(left_len > 0) {
    	for(int i = count; i < sss.size(); i ++) {
    		temp <<= 1;
    		if(sss[i] == '1') temp |= 1;
		}
		for(int i = 0; i < 8 - left_len; i ++) {
			temp <<= 1;
			temp |= 1;
		}
	}
    result.push_back(temp);

    return result;
}

vector<unsigned char> huffman_decompress(vector<unsigned char> binCode){
    huf_tree test;
	vector<unsigned char> result = {};
    

    int file_len = binCode.size();
    unsigned char s;
    char st[8];
    unsigned char tmp;
	char arr[file_len];

    
	for (int i=0;i<file_len;i++)
        arr[i] = binCode[i];
  
    char storage[10 * file_len];
    int count = 0, cur = 0;
    for(int j = 0; j < file_len; j ++) {
		unsigned char c = arr[j];

		char a[8];
		count += 8;
 		for(int i = 1; i <= 8; i ++) {
			storage[count - i] = (c % 2) + '0';
			c /= 2;
		}
	}
	for(int i = 0; i < 8; i ++) st[i] = storage[cur ++];
	
    s = chTobit(st);

	//cout << "Decompressed byte "<<int(s) << endl;
    result.push_back(s);
    test.addnode(s);// 
    while(1) {
        node *p = new node;
        p = test.root;
        unsigned char c;
        while(1) {
			c = storage[cur ++];
            if(c == '0') {
                p = p->left;
            }
            else {
                p = p->right;
            }
            if(p->uch != 250 && p->num != 0) {
                s = p->uch;
                break;
            }
            if(p->uch == 250 && p->left == NULL && p->right == NULL) {
                s = 250;
                break;
            }
        }
        if(s != 250) {
            result.push_back(s);
            test.addnode(s);
        }
        else {
			for(int i = 0; i < 8; i ++) st[i] = storage[cur ++];
            s = chTobit(st);
            if(s == 255) break;
            result.push_back(s);

            test.addnode(s);
        }
    }
    
    return result;
}

vector<vector<unsigned char>> I_frame_encoder (vector<vector<unsigned char>> image,vector<vector<double>> quanti_matrix){
    int num_block_height = int(ceil(double(image.size())/block_size));
    int num_block_width = int(ceil(double(image[0].size())/block_size));
    int padded_height = num_block_height * block_size;
    int padded_width = num_block_width * block_size;

    // add padding to the image so that the size is multiples of block_size
    vector<vector<double>> padded_image(padded_height, vector<double>(padded_width, 0));
    //vector<vector<double>> padded_image(padded_height, vector<double>(padded_width, 0));
    for(int y = 0; y < image.size(); y++){
        for (int x = 0; x < image[0].size(); x++){
            padded_image[y][x] =  image[y][x] - 128;
        }
    }
    for (int y = 0; y < image.size();y++){
        for (int x = image[0].size(); x < padded_image[0].size();x++){
            padded_image[y][x] = image[y][image[0].size()-1] - 128;
        }
    }
    
    for (int x = 0; x < padded_image[0].size();x++){
        for (int y = image.size(); y < padded_image.size();y++){
            padded_image[y][x] = image[image.size()-1][x] - 128;
        }
    }
    
    auto block = create_2d_vector<double>(int(block_size),int(block_size));
    


    for (int col_idx = 0; col_idx < num_block_height ; col_idx++){
        int col_idx_start = col_idx * block_size;
        
        for (int row_idx = 0; row_idx < num_block_width  ; row_idx++){
            int row_idx_start = row_idx * block_size;

            // extract the contents in the block
            for (int i = 0;i<block_size;i++){
                for (int j =0;j<block_size;j++){
                    block[i][j] = padded_image[col_idx_start + i][row_idx_start + j];
                }
            }
            
     
            // DCT transform
            block = dctTransform(block);


  


            // Quantization transform
            for (size_t i = 0; i < block_size; i++) {
                for (size_t j = 0; j < block_size; j++) {
                    block[i][j] = round(block[i][j] / quanti_matrix[i][j]);
                }
            }


  
            
            vector<vector<unsigned char>> block_downsampled(8, vector<unsigned char>(8, 128));
            for (auto i: block_cords){
                int x_index = i[0];
                int y_index = i[1];
                block_downsampled[x_index][y_index] = block[x_index][y_index] + 128;
            }         


            auto block_flattened = zig_zag_flatten(block_downsampled);
            auto block_delta = delta_compress(block_flattened);
            auto block_huffman = huffman_compress(block_delta);*/
            

            unsigned char DC_coeff = block_downsampled[0][0];
            auto t1 = RLE_without_first_element(block_downsampled);
            //auto t2 = huffman_compress(t1);
            auto t2 = t1;




            if (write_data){
                output_stream.push_bytes(DC_coeff);
                output_stream.push_byte(t2.size());
                for (auto byte:t2){
                    output_stream.push_byte(byte);
                }
            }
            
            
            //cout << block_flattened.size()<<endl;
            //cout << block_delta.size()<<endl;
            //cout << block_huffman.size()<<endl;
            
           
            
            
            
            // decompressed frame 

          

            for (size_t i = 0; i < block_size; i++) {
                for (size_t j = 0; j < block_size; j++) {
                    block[i][j] = round(block[i][j] * quanti_matrix[i][j]);
                }
            }

    
            block = inverse_dctTransform(block);


    
            
            for (int i =0;i<8;i++){
                for (int j=0;j<8;j++){
                    block[i][j] = round_and_clamp_to_char(block[i][j] + 128); 
                }
            }

            // copy blocks back to image 
            for (int i = 0;i<block_size;i++){
                for (int j=0;j<block_size;j++){
                    padded_image[col_idx_start+i][row_idx_start+j] = block[i][j];
                }
            }     
        }
    }
    for (int i = 0; i<image.size();i++){
        for (int j =0;j<image[0].size();j++){
            image[i][j] = padded_image[i][j];
        }
    }
    return image;
}

int motion_vector_evaluate(vector<vector<unsigned char>>image,vector<vector<unsigned char>>reference,int division_factor){
    int sum = 0;
    for (int i=0;i<macro_size;i++){
        for (int j=0;j<macro_size;j++){
            sum = sum + abs(image[i][j] - reference[i][j]);
        }
    };
    return sum/division_factor;
}

// padding the image to multiple of 16 and shift range from 0-255 to -128 to 127
vector<vector<unsigned char>> image_padder(vector<vector<unsigned char>> image,int size){
    int num_block_height = int(ceil(double(image.size())/size));
    int num_block_width = int(ceil(double(image[0].size())/size));
    int padded_height = num_block_height * size;
    int padded_width = num_block_width * size;
    vector<vector<unsigned char>> padded_image(padded_height, vector<unsigned char>(padded_width, 0));
    for(int y = 0; y < image.size(); y++){
        for (int x = 0; x < image[0].size(); x++){
            padded_image[y][x] =  image[y][x];
        }
    }
   return padded_image;
} 

// find the best matching of current macro block and the reference frame (requires height and width to be multiple of 16)
vector<int> motion_vector_finder(vector<vector<unsigned char>> current_macro, vector<vector<unsigned char>> reference_frame, int current_idx_y, int current_idx_x){
    int height_start = max(0, current_idx_y - search_range);
    int height_end = min(int(reference_frame.size() -macro_size), current_idx_y + search_range);
    int width_start = max(0, current_idx_x - search_range);
    int width_end = min(int(reference_frame[0].size()  - macro_size), current_idx_x + search_range);
    vector<vector<unsigned char>> reference_block(macro_size,vector<unsigned char>(macro_size,0));
    int score;
    int height_offset = 0;
    int width_offset = 0;
    int best_score = 10000;
    //cout << height_start << " - " << height_end << " || " << width_start << " - " << width_end<<endl;
    for (int i=height_start;i<height_end;i=i+1){
        for (int j=width_start;j<width_end;j=j+1){ 
            // Extrac the reference block
            for (int m=0;m<macro_size;m++){
                for (int n=0;n<macro_size;n++){
                    //cout << i+m << " - " << j + n << endl;
                    reference_block[m][n] = reference_frame[i+m][j+n];
                }
            }
            score =  motion_vector_evaluate(current_macro,reference_block,macro_size^2);
            if (score < best_score){
                height_offset = i - current_idx_y;
                width_offset = j - current_idx_x;
                best_score = score;
            }
        }
    }
    return {height_offset,width_offset};
}

vector<vector<int>> motion_vectors_finder(vector<vector<unsigned char>> current_frame,vector<vector<unsigned char>> reference_frame){
    vector<vector<int>> motion_vectors = {};
    
    vector<vector<unsigned char>> current_block = create_2d_vector<unsigned char>(macro_size,macro_size);
    
    vector<vector<unsigned char>> current_frame_padded = image_padder(current_frame,16);
    vector<vector<unsigned char>> reference_frame_padded = image_padder(reference_frame,16);

    int num_block_height = int(ceil(double(current_frame_padded.size())/macro_size));
    int num_block_width = int(ceil(double(current_frame_padded[0].size())/macro_size));
 

    for (int i = 0;i < num_block_height ; i++){
        for (int j = 0; j < num_block_width ; j++){


            // index of the block 
            int macro_block_height = i * macro_size;
            int macro_block_width  = j * macro_size;


            // extract blocks from the current frame
            for (int m=0;m<macro_size;m++){
                for (int n=0;n<macro_size;n++){
                    //cout <<  macro_block_height + m << " - " << macro_block_width + n << endl;
                    current_block[m][n] = current_frame_padded[macro_block_height + m][macro_block_width + n];
                }
            }
            auto motion_vector = motion_vector_finder(current_block,reference_frame_padded,macro_block_height,macro_block_width);
            
     
            motion_vectors.push_back(motion_vector);
        }
    }
    
    

    return motion_vectors;
}

vector<vector<unsigned char>> frame_difference(vector<vector<unsigned char>> image, vector<vector<unsigned char>> ref, vector<vector<int>> MVs){
    auto image_padded = image_padder(image,16);
    auto ref_padded = image_padder(ref,16);
    vector<vector<unsigned char>> diff(image_padded.size(), vector<unsigned char>(image_padded[0].size(), 0));
    int num_block_height = int(ceil(double(image_padded.size())/macro_size));
    int num_block_width = int(ceil(double(image_padded[0].size())/macro_size));
    int block_counter = 0;
    for (int i = 0;i < num_block_height ; i++){
        for (int j = 0; j < num_block_width ; j++){
            int height_offset = MVs[block_counter][0];
            int width_offset = MVs[block_counter][1];

            // index of the block 
            int macro_block_height = i * macro_size;
            int macro_block_width  = j * macro_size;

            for (int m=0;m<macro_size;m++){
                for (int n=0;n<macro_size;n++){
                    int idx_height = macro_block_height+m;
                    int idx_width = macro_block_width+n;
                    //cout <<  idx_height << " - " << height_offset << " | " << idx_width << " - " << width_offset << endl; 
                    diff[idx_height][idx_width] = (image[idx_height][idx_width] - ref[idx_height+height_offset][idx_width+width_offset])/2 + 128;
                    //cout << "p 2" << endl;
                }
            }
            block_counter ++;
        }
    }
    return diff;
}

vector<vector<unsigned char>> frame_from_difference(vector<vector<unsigned char>> diff, vector<vector<unsigned char>> ref, vector<vector<int>> MVs ){
    int num_block_height = int(ceil(double(diff.size())/macro_size));
    int num_block_width = int(ceil(double(diff[0].size())/macro_size));
    //cout << num_block_height << endl;
    //cout << num_block_width << endl;
    auto result = create_2d_vector<unsigned char>(diff.size(),diff[0].size());
    //cout << MVs.size()<<endl;
    /*
    for (int i=0;i<MVs.size();i++){
        auto pair = MVs[i];
        cout << i << " < ";
        cout << pair[0] << " || " << pair[1] << endl;
    }cout << endl;*/
    
    for (int i=0;i<diff.size();i++){
        for (int j=0;j<diff[0].size();j++){
            int block_idx_height = i/macro_size;
            int block_idx_width = j/macro_size;
            //cout << i << " - " << j << endl;
            //cout << block_idx_height << " - "<<block_idx_width<<endl;
            int order = block_idx_height*num_block_width+block_idx_width;
            //out << block_idx_height << " " << num_block_width << " " << block_idx_width<<endl;
            //cout << "1" << endl;
            //cout << order << endl<<endl;
            int height_offset = MVs[order][0];
            int width_offset = MVs[order][1];
            //cout << "order " << order << " ";
            //cout << height_offset << " "<<width_offset << endl;

            //cout << "3" << endl;
            //cout <<  i << " - " << height_offset << " | " << j << " - " << width_offset << endl; 
            //cout << int(diff[i][j])<<endl;
            result[i][j] = round_and_clamp_to_char((int(diff[i][j])-128) * 2 + ref[i+height_offset][j+width_offset]);
            //cout << "4" << endl;

        }
    }
    return result;
}

vector<vector<unsigned char>> P_frame_encoder(vector<vector<unsigned char>> image, vector<vector<unsigned char>> ref, vector<vector<double>> quanti_matrix){

    auto MVs = motion_vectors_finder(image,ref);
   
    if (write_data){
        for (int i =0;i<MVs.size();i++){  
            output_stream.push_byte(MVs[i][0] + 128);
            output_stream.push_byte(MVs[i][1] + 128);  
        }      
    }
    

    
    
    
    vector<vector<unsigned char>> diff = frame_difference(image,ref,MVs);
    vector<vector<unsigned char>> diff_decompressed = I_frame_encoder(diff,quanti_matrix);
    auto res = frame_from_difference(diff_decompressed,ref,MVs);


    return res;
}

int main(int argc, char** argv){
    
    
    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <width> <height> <low/medium/high>" << std::endl;
        return 1;
    }
    u32 width = std::stoi(argv[1]);
    u32 height = std::stoi(argv[2]);
    std::string quality{argv[3]};

    YUVStreamReader reader {std::cin, width, height};
    OutputBitStream output_stream {std::cout};
    vector<vector<unsigned char>> frame_Y,frame_Cb,frame_Cr;
    output_stream.push_u32(height);
    output_stream.push_u32(width);

    double quanti_multi;
    if (quality == "low"){
        quanti_multi = 10;
        output_stream.push_bytes(1);

    }else if (quality == "medium"){
        quanti_multi = 5;
        output_stream.push_bytes(100);

    }else if (quality == "high"){
        quanti_multi = 1;
        output_stream.push_bytes(200);

    }else if (quality == "full"){
        quanti_multi = 0.5;
        output_stream.push_bytes(255);
    }else{
        quanti_multi = 1;
    }
 
    for (int i=0;i<block_size;i++){
        for (int j=0;j<block_size;j++){
            if (i + j > DCT_Tune){
                quantization_matrix_lumin[i][j]*= quanti_multi;
                quantization_matrix_chrom[i][j]*= quanti_multi;

            }else{
                continue;
            }
        }
    }

    block_cords = zig_zag_coords(8,8,64);
    while (reader.read_next_frame()){
        output_stream.push_byte(1); //Use a one byte flag to indicate whether there is a frame here
        YUVFrame420& frame = reader.frame();


        unsigned int scaled_height = (height)/2;
        unsigned int scaled_width = (width)/2;
        auto Y = create_2d_vector<unsigned char>(height,width);
        auto Cb = create_2d_vector<unsigned char>(scaled_height,scaled_width);
        auto Cr = create_2d_vector<unsigned char>(scaled_height,scaled_width);
        
        for (u32 y = 0; y < height; y++)
            for (u32 x = 0; x < width; x++)
                Y[y][x] = frame.Y(x,y);
        for (u32 y = 0; y < height/2; y++)
            for (u32 x = 0; x < width/2; x++)
                Cb[y][x] = frame.Cb(x,y);
        for (u32 y = 0; y < height/2; y++)
            for (u32 x = 0; x < width/2; x++)
                Cr[y][x] = frame.Cr(x,y);
        
        if (counter == 0){
            frame_Y = I_frame_encoder(Y,quantization_matrix_lumin);
            frame_Cb = I_frame_encoder(Cb,quantization_matrix_chrom);
            frame_Cr = I_frame_encoder(Cr,quantization_matrix_chrom);    
        
            Y_buffer = (frame_Y);
            Cb_buffer = (frame_Cb);
            Cr_buffer = (frame_Cr);
            


        }else{
    
            
            frame_Y  = P_frame_encoder(Y,Y_buffer,quantization_matrix_lumin);
            frame_Cb = P_frame_encoder(Cb,Cb_buffer,quantization_matrix_chrom);
            frame_Cr = P_frame_encoder(Cr,Cr_buffer,quantization_matrix_chrom);
            Y_buffer = (frame_Y);
            Cb_buffer = (frame_Cb);
            Cr_buffer = (frame_Cr);
        }  
      
        
     
        
    
        counter++;
        total_frames ++ ;
        //cout << "Frame Counter "<<total_frames << endl;
        if (counter == 60){
            counter = 0;  
            
        }
    }
    
    output_stream.push_byte(0); //Flag to indicate end of data
    output_stream.flush_to_byte();

    return 0;
}
