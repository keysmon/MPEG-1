----Overview----

MPEG-1 implementation in C++. It uses I frames and P frames. I frames are applied to every 15 frames and P frames are used for the other 14 frames. 
Transformation: RLE, Delta and adaptive Huffman compression. Delta compression was implemented but not used due to less compression ratio in testing. The quality setting is achieved by adjusting the quantization matrix and the percentage of DCT coefficients to keep. 

----Feature Implemented----

BASIC
1. DOCUMENTATION
2. P-frames implementation
3. motion compensation
4. 12 compression ratio

ADVANCED
1. Adaptive Huffman coding

----Architecture----
I Frames:
1. Divide each frame into 8 by 8 blocks. 
2. Shift pixel range down by 128
3. DCT and quantization
4. Shift range up by 128
5. RLE and adaptive Huffman coding

P Frames:
1. Divide each frame into 16 by 16 blocks
2. Find motion vectors using a local search
3. Find the differences between the current frame and the previous frame with motion vectors applied
4. Apply I-frame transformation to store differences.




----Bibliography----

1. huffman_tree.hpp are sourced from (https://download.csdn.net/download/qq_40365272/12104576?ops_request_misc=%257B%2522request%255Fid%2522%253A%2522169259773016800185833201%2522%252C%2522scm%2522%253A%252220140713.130102334.pc%255Fall.%2522%257D&request_id=169259773016800185833201&biz_id=1&utm_medium=distribute.pc_search_result.none-task-download-2~all~first_rank_ecpm_v1~rank_v31_ecpm-3-12104576-null-null.142^v93^insert_down28v1&utm_term=%E8%87%AA%E9%80%82%E5%BA%94%E5%93%88%E5%A4%AB%E6%9B%BCc%2B%2B&spm=1018.2226.3001.4187.4). It contains the general structure of an adaptive huffman tree. I implemented encoded and decoding functions using this tree. 










