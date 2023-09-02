/* uvid_decompress.cpp
   CSC 485B/578B - Data Compression - Summer 2023

   Starter code for Assignment 4
   
   This placeholder code reads the (basically uncompressed) data produced by
   the uvid_compress starter code and outputs it in the uncompressed 
   YCbCr (YUV) format used for the sample video input files. To play the 
   the decompressed data stream directly, you can pipe the output of this
   program to the ffplay program, with a command like 

     ffplay -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 - 2>/dev/null
   (where the resolution is explicitly given as an argument to ffplay).

   B. Bird - 2023-07-08
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include <tuple>
#include "input_stream.hpp"
#include "yuv_stream.hpp"
#include <math.h> 
#include "output_stream.hpp"
#include "uvg_common.hpp"
#include "huffman_tree.hpp"
#include "cmath"


#define pi 3.1415926
const int block_size  = 8;
const int macro_size = 16;
const int DCT_Tune = 3;
u32 height;
u32 width;
int counter = 0;
int total_frames = 0;
int num_block_height_Y,num_block_width_Y,num_block_height_Cbr,num_block_width_Cbr;
int padded_height_Y,padded_width_Y,padded_height_Cbr,padded_width_Cbr;
vector<vector<unsigned char>> Y_prev;
vector<vector<unsigned char>> Cb_prev;
vector<vector<unsigned char>> Cr_prev;

void num_of_block_calculator(int height,int width,int &num_block_height_Y,int &num_block_width_Y,int &num_block_height_Cbr,int &num_block_width_Cbr){
    num_block_height_Y = int(ceil(double(height)/block_size));
    num_block_width_Y = int(ceil(double(width)/block_size));
    padded_height_Y = num_block_height_Y * block_size;
    padded_width_Y = num_block_width_Y * block_size;

    int Cbr_height = (height+1)/2;
    int Cbr_width = (width+1)/2;
    
    num_block_height_Cbr = int(ceil(double(Cbr_height)/block_size));
    num_block_width_Cbr = int(ceil(double(Cbr_width)/block_size));
    padded_height_Cbr = num_block_height_Cbr * block_size;
    padded_width_Cbr = num_block_width_Cbr * block_size;
}

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

vector<vector<double>> inverse_dctTransform(vector<vector<double>> matrix)
{   
    
    int u,v,x,y;
   
    // dct will store the discrete cosine transform
    vector<vector<double>> dct = create_2d_vector<double>(block_size,block_size);
    double au, av, dct1, sum;
    for (x = 0;  x < block_size; x++) {
        for (y = 0; y < block_size; y++) {
            // ci and cj depends on frequency as well as
            // number of row and columns of specified matrix
            
            // sum will temporarily store the sum of
            // cosine signals
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

vector<unsigned char> huffman_decompress(vector<unsigned char> binCode){
    huf_tree test;
	vector<unsigned char> result = {};
    

    int file_len = binCode.size();
    unsigned char s;
    char st[8];
    unsigned char tmp;
	char arr[file_len];

    
    //cout<<"file_len:"<<file_len<<endl;
    //cout<<"decompressing........."<<endl;
	
	for (int i=0;i<file_len;i++)
        arr[i] = binCode[i];
    //for (auto byte:arr)
        //cout << int(byte) << endl;
	

	
	
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
			//cout << "Decompressed byte "<<int(s) << endl;
            result.push_back(s);

            test.addnode(s);
        }
        else {
			for(int i = 0; i < 8; i ++) st[i] = storage[cur ++];
            s = chTobit(st);
            if(s == 255) break;
			//cout << "Decompressed byte "<<int(s) << endl;
            result.push_back(s);

            test.addnode(s);
        }
    }
    /*
    for (auto i:result){
        cout << i;
    }cout << endl;*/
    return result;
}

vector<unsigned char> inverse_RLE(vector<unsigned char> rle){
    vector<unsigned char> res = {};
    for (int i=0;i<rle.size();i++){
        for (int j=0;j<rle[i+1];j++){
            res.push_back(rle[i]);
        }i = i + 1;
    }return res;
}

vector<vector<unsigned char>> I_frame_decoder(vector<unsigned char> DC_coeffs,vector<vector<unsigned char>> seqeuence, vector<vector<double>> quanti_matrix, int original){
    int num_block_height,num_block_width;
    //cout << "p1" << endl;
    if (original == 1){
        num_block_height = num_block_height_Y;
        num_block_width = num_block_width_Y; 
    }else{
        num_block_height = num_block_height_Cbr;
        num_block_width = num_block_width_Cbr;
    }
    int padded_height = num_block_height * block_size;
    int padded_width = num_block_width * block_size;
    vector<vector<unsigned char>> padded_image(padded_height, vector<unsigned char>(padded_width, 0));
    //cout << "p2" << endl;
    for (int s = 0;s<seqeuence.size();s++){

        vector<vector<double>> block(block_size, vector<double>(block_size, 0));
        int block_idx_height = s/num_block_width;
        int block_idx_width = s%num_block_width;

        auto seq = seqeuence[s];
        //auto huffman_decoded = huffman_decompress(seq);
        auto huffman_decoded = seq;



        auto rle_decoded = inverse_RLE(huffman_decoded);
        rle_decoded.insert(rle_decoded.begin(),DC_coeffs[s]);
        
        /*cout << "-- " << int(DC_coeffs[s]) << " --" << endl;
        for (auto i:seq){
            cout << int(i) << " ";
        }cout << endl;
        for (auto i:huffman_decoded){
            cout << int(i) << " ";
        }cout << endl<<endl;*/
        
        

        
        for (int i=0;i<64 ;i++){
            block[zig_zag_route[i][0]][zig_zag_route[i][1]] = double(rle_decoded[i]);
        }
       

        
        /*for (auto i:seq){
            cout << int(i) << " " ;
        }cout << endl<< endl;
        
        for (auto i:huffman_decoded){
            cout << int(i) << " " ;
        }cout << endl<< endl;
        
        for (auto i:rle_decoded){
            cout << int(i) << " " ;
        }cout << endl<< endl;
        
        cout << endl << " ------ " << endl;
        for (int i=0;i<block_size;i++){
            for (int j=0;j<block_size;j++){
                cout << block[i][j] << " ";
            }cout << endl;
        }cout << endl<<endl;*/
        
        

        // decompressed frame 
        for (size_t i = 0; i < block_size; i++) {
            for (size_t j = 0; j < block_size; j++) {
                block[i][j] = round((block[i][j]-128) * quanti_matrix[i][j]);
            }
        }
        
        /*cout << "Dequantized " << endl;
        for (int i=0;i<block_size;i++){
            for (int j=0;j<block_size;j++){
                cout << block[i][j] << " ";
            }cout << endl;
        }cout << endl<<endl;*/


        block = inverse_dctTransform(block);
        

        /*
        cout << "Reverse DCT" << endl;
        for (int i=0;i<block_size;i++){
            for (int j=0;j<block_size;j++){
                cout << block[i][j] << " ";
            }cout << endl;
        }cout << endl<<endl;*/
        
        for (int i =0;i<8;i++){
            for (int j=0;j<8;j++){
                block[i][j] = round_and_clamp_to_char(block[i][j] + 128); 
            }
        }
        
        /*
        cout << "Decompressed result" << endl;
        for (int i=0;i<block_size;i++){
            for (int j=0;j<block_size;j++){
                cout << block[i][j] << " ";
            }cout << endl;
        }cout << endl<<endl;*/
        
        for (int i=0;i<block_size;i++){
            for (int j=0;j<block_size;j++){
                padded_image[block_idx_height*block_size+i][block_idx_width*block_size+j] = block[i][j];
            }
        }
        //cout << "Fill it back - "<< s << endl;
    }
    //cout << "All blocks filled back" << endl;
    
    /*
    for (int i=0;i<10;i++){
        for (int j=0;j<10;j++){
            cout << int(padded_image[i][j]) << " ";
        }cout << endl;
    }cout << endl<<endl;*/
    
    return padded_image;
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

            int order = block_idx_height*num_block_width+block_idx_width;           
             //cout << "1" << endl;

            int height_offset = MVs[order][0];
            int width_offset = MVs[order][1];
            //cout << "order " << order << " ";
            //cout << height_offset << " "<<width_offset << endl;

            //cout << "3" << endl;
            //cout <<  i << " - " << height_offset << " | " << j << " - " << width_offset << endl; 

            result[i][j] = round_and_clamp_to_char((int(diff[i][j])-128) * 2 + ref[i+height_offset][j+width_offset]);
            
            //cout << "4" << endl;

        }
    }





    return result;
}

vector<vector<unsigned char>> P_frame_decoder(vector<vector<unsigned char>> MVs, vector<unsigned char> DC_coeffs,vector<vector<unsigned char>> seqeuence, vector<vector<unsigned char>> ref,vector<vector<double>> quanti_matrix,int original){

    
    
    vector<vector<int>> MVs_shifted(MVs.size(), vector<int>(2, 0));
    for (int i = 0;i<MVs.size();i++){
        MVs_shifted[i][0] = int(MVs[i][0]) - 128;
        MVs_shifted[i][1] = int(MVs[i][1]) - 128;
        //cout << int(MVs[i][0])  << " || " << int(MVs[i][1])  << endl;
    }

    /*
    cout << endl;
    for (int i = 0;i<DC_coeffs.size();i++){
        
        cout << int(DC_coeffs[i])  << " || " << int(DC_coeffs[i])  << endl;
    }*/



    auto diff = I_frame_decoder(DC_coeffs,seqeuence,quanti_matrix,original);
    auto res = frame_from_difference(diff,ref,MVs_shifted);

    /*
    for (int i=0;i<ref.size();i++){
        for (int j=0;j<ref[0].size();j++){
            res[i][j] = round_and_clamp_to_char((diff[i][j] - 128) * 2 + ref[i][j]);
        }
    }*/
    



    return res;
    
}


int main(int argc, char** argv){

    //Note: This program must not take any command line arguments. (Anything
    //      it needs to know about the data must be encoded into the bitstream)
    
    
    InputBitStream input_stream {std::cin};

    height =  input_stream.read_u32();
    width = input_stream.read_u32();

    int num_macro_height_Y = int(ceil(double(height)/macro_size));
    int num_macro_width_Y = int(ceil(double(width)/macro_size));
    int num_of_MVs_Y = num_macro_height_Y * num_macro_width_Y ;

    int num_macro_height_Cbr = int(ceil(double(height/2)/macro_size));
    int num_macro_width_Cbr = int(ceil(double(width/2)/macro_size)); 
    int num_of_MVs_Cbr = num_macro_height_Cbr * num_macro_width_Cbr;
    

    num_of_block_calculator(height, width, num_block_height_Y,num_block_width_Y,num_block_height_Cbr,num_block_width_Cbr);

    unsigned char quality  = input_stream.read_byte();
    double quanti_multi = 0.5;
    if (quality == 1){
        quanti_multi = 10;
    }else if (quality == 100){
        quanti_multi = 5;
    }else if (quality == 200){
        quanti_multi = 1;
    }else if (quality == 255){
        quanti_multi = 0.5;
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


    YUVStreamWriter writer {std::cout, width, height};
    vector<vector<unsigned char>> buffer =  {};
    vector<unsigned char> DC_coeffs = {};
    vector<unsigned char> byte;
    vector<vector<unsigned char>> MVs;
    while (input_stream.read_byte()){
        YUVFrame420& frame = writer.frame();
        vector<vector<unsigned char>> Y_frame,Cb_frame,Cr_frame;
        
        // Identity Frames
        if (counter == 0){

            //read the bytes
            for (int a = 0;a<num_block_height_Y*num_block_width_Y;a++){
                unsigned char DC_coeff = input_stream.read_byte();
                DC_coeffs.push_back(DC_coeff);
                int num_of_byte = input_stream.read_byte();
                for (int i=0;i<num_of_byte;i++){
                    byte.push_back(input_stream.read_byte());
                }
                buffer.push_back(byte);
                byte.clear();
            }
            Y_frame = I_frame_decoder(DC_coeffs,buffer,quantization_matrix_lumin,1);
            buffer.clear();//clear buffer
            DC_coeffs.clear();






            //read the bytes
            for (int a = 0;a<num_block_height_Cbr*num_block_width_Cbr;a++){
                unsigned char DC_coeff = input_stream.read_byte();
                DC_coeffs.push_back(DC_coeff);
                int num_of_byte = input_stream.read_byte();
                for (int i=0;i<num_of_byte;i++){
                    byte.push_back(input_stream.read_byte());
                }
                buffer.push_back(byte);
                byte.clear();
            }
            Cb_frame = I_frame_decoder(DC_coeffs,buffer,quantization_matrix_chrom,0);
            buffer.clear();//clear buffer
            DC_coeffs.clear();


            
            //read the bytes
            for (int a = 0;a<num_block_height_Cbr*num_block_width_Cbr;a++){
                unsigned char DC_coeff = input_stream.read_byte();
                DC_coeffs.push_back(DC_coeff);
                int num_of_byte = input_stream.read_byte();
                
                for (int i=0;i<num_of_byte;i++){
                    byte.push_back(input_stream.read_byte());
                }
                buffer.push_back(byte);
                byte.clear();
            }
            Cr_frame = I_frame_decoder(DC_coeffs,buffer,quantization_matrix_chrom,0);
            buffer.clear();//clear buffer
            DC_coeffs.clear();

            Y_prev = (Y_frame);
            Cb_prev = (Cb_frame);
            Cr_prev = (Cr_frame);
           

        // Prediction_frame
        }else{
            //read the bytes
            for (int i = 0 ; i < num_of_MVs_Y;i = i+1){
                unsigned char temp1 = input_stream.read_byte();
                unsigned char temp2 = input_stream.read_byte();
                MVs.push_back({temp1,temp2});
            }
            //cout << MVs.size()<<endl;
            for (int a = 0;a<num_block_height_Y*num_block_width_Y;a++){
                unsigned char DC_coeff = input_stream.read_byte();
                DC_coeffs.push_back(DC_coeff);
                int num_of_byte = input_stream.read_byte();
                for (int i=0;i<num_of_byte;i++){
                    byte.push_back(input_stream.read_byte());
                }
                buffer.push_back(byte);
                byte.clear();
            }
            Y_frame = P_frame_decoder(MVs,DC_coeffs,buffer,Y_prev,quantization_matrix_lumin,1);
            buffer.clear();//clear buffer
            DC_coeffs.clear();
            MVs.clear();





            for (int i = 0 ; i < num_of_MVs_Cbr;i = i+1){
                unsigned char temp1 = input_stream.read_byte();
                unsigned char temp2 = input_stream.read_byte();
                MVs.push_back({temp1,temp2});
            }
            for (int a = 0;a<num_block_height_Cbr*num_block_width_Cbr;a++){
                unsigned char DC_coeff = input_stream.read_byte();
                DC_coeffs.push_back(DC_coeff);
                int num_of_byte = input_stream.read_byte();
                for (int i=0;i<num_of_byte;i++){
                    byte.push_back(input_stream.read_byte());
                }
                buffer.push_back(byte);
                byte.clear();
            }
            Cb_frame = P_frame_decoder(MVs,DC_coeffs,buffer,Cb_prev,quantization_matrix_chrom,0);
            buffer.clear();//clear buffer
            DC_coeffs.clear();
            MVs.clear();







            for (int i = 0 ; i < num_of_MVs_Cbr;i = i+1){
                unsigned char temp1 = input_stream.read_byte();
                unsigned char temp2 = input_stream.read_byte();
                MVs.push_back({temp1,temp2});
            }
            for (int a = 0;a<num_block_height_Cbr*num_block_width_Cbr;a++){
                unsigned char DC_coeff = input_stream.read_byte();
                DC_coeffs.push_back(DC_coeff);
                int num_of_byte = input_stream.read_byte();
                for (int i=0;i<num_of_byte;i++){
                    byte.push_back(input_stream.read_byte());
                }
                buffer.push_back(byte);
                byte.clear();
            }
            Cr_frame = P_frame_decoder(MVs,DC_coeffs,buffer,Cr_prev,quantization_matrix_chrom,0);
            buffer.clear();//clear buffer
            DC_coeffs.clear();
            MVs.clear();

            Y_prev = (Y_frame);
            Cb_prev = (Cb_frame);
            Cr_prev = (Cr_frame);
            
        }
        /*
        //cout << "---- frame  "<< total_frames<<" ----" << endl;
        for (int i = 0;i<10;i++){
            for (int j=0;j<10;j++){
                cout << int(Y_frame[i][j]) << " ";
            }cout << endl;
        }cout << endl;
        for (int i = 0;i<10;i++){
            for (int j=0;j<10;j++){
                cout << int(Cb_frame[i][j]) << " ";
            }cout << endl;
        }cout << endl;
        for (int i = 0;i<10;i++){
            for (int j=0;j<10;j++){
                cout << int(Cr_frame[i][j]) << " ";
            }cout << endl;
        }cout << endl<<endl;*/

        counter ++;
        total_frames ++;
        //cout << total_frames << endl;
        if (counter == 60){
            counter = 0;
            
        }
        //cout << total_frames << endl;
        for (u32 y = 0; y < height; y++)
            for (u32 x = 0; x < width; x++)
                frame.Y(x,y) = Y_frame[y][x];
        for (u32 y = 0; y < height/2; y++)
            for (u32 x = 0; x < width/2; x++)
                frame.Cb(x,y) = Cb_frame[y][x];
        for (u32 y = 0; y < height/2; y++)
            for (u32 x = 0; x < width/2; x++)
                frame.Cr(x,y) = Cr_frame[y][x];
        
        writer.write_frame();
        
    }




    return 0;
}