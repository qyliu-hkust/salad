#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <bitset>
#include <fstream>
#include <sstream>

// 哈夫曼树节点定义
struct HuffmanNode {
    int value;  // 存储符号的值
    int freq;   // 符号出现的频率
    HuffmanNode* left;
    HuffmanNode* right;

    HuffmanNode(int v, int f) : value(v), freq(f), left(nullptr), right(nullptr) {}

    // 用于优先队列的比较器
    bool operator>(const HuffmanNode& other) const {
        return freq > other.freq;  // 按频率从小到大排序
    }
};

// 手动实现的最小堆优先队列
class MinHeap {
public:
    // 向堆中插入元素
    void push(HuffmanNode* node) {
        heap.push_back(node);  // 插入到堆的末尾
        bubbleUp(heap.size() - 1);  // 上浮到正确的位置
    }

    // 移除并返回堆顶元素
    HuffmanNode* pop() {
        if (heap.empty()) return nullptr;

        // 将堆顶元素与最后一个元素交换
        std::swap(heap[0], heap[heap.size() - 1]);
        HuffmanNode* top = heap.back();
        heap.pop_back();

        // 下沉根元素到正确位置
        if (!heap.empty()) {
            bubbleDown(0);
        }

        return top;
    }

    HuffmanNode* top() const {
        return heap.empty() ? nullptr : heap[0];
    }

    // 返回堆的副本，用于打印堆的内容
    std::vector<HuffmanNode*> getHeap() const {
        std::vector<HuffmanNode*> heap_copy(heap);
        return heap_copy;
    }

    // 返回堆的大小
    size_t size() const {
        return heap.size();
    }

private:
    std::vector<HuffmanNode*> heap;

    // 上浮操作：将当前元素上浮到正确的位置
    void bubbleUp(size_t index) {
        while (index > 0) {
            size_t parent = (index - 1) / 2;
            if (*heap[index] > *heap[parent]) break;  // 如果当前元素大于父节点，则不需要上浮

            // 否则，交换当前元素和父节点
            std::swap(heap[index], heap[parent]);
            index = parent;
        }
    }

    // 下沉操作：将当前元素下沉到正确的位置
    void bubbleDown(size_t index) {
        size_t leftChild, rightChild, smallest;

        while (index < heap.size()) {
            leftChild = 2 * index + 1;
            rightChild = 2 * index + 2;
            smallest = index;

            // 找到最小的子节点
            if (leftChild < heap.size() && heap[leftChild]->freq < heap[smallest]->freq) {
                smallest = leftChild;
            }
            if (rightChild < heap.size() && heap[rightChild]->freq < heap[smallest]->freq) {
                smallest = rightChild;
            }

            // 如果当前节点已经是最小的，则停止
            if (smallest == index) break;

            // 否则交换当前节点和最小子节点
            std::swap(heap[index], heap[smallest]);
            index = smallest;
        }
    }
};

// 用于构建哈夫曼树的函数
HuffmanNode* buildHuffmanTree(const std::unordered_map<int, int>& freqMap) {
    // 创建手动实现的最小堆
    MinHeap minHeap;

    // 创建哈夫曼树的初始节点并加入堆中
    std::cout << "Initial Nodes (based on frequencies):\n";
    for (const auto& entry : freqMap) {
        HuffmanNode* node = new HuffmanNode(entry.first, entry.second);
        minHeap.push(node);
        std::cout << "Node: Value = " << node->value << ", Frequency = " << node->freq << std::endl;
    }

    // 构建哈夫曼树
    std::cout << "\nBuilding Huffman Tree:\n";
    while (minHeap.size() > 1) {
        // 每次从堆中取出两个最小的节点
        HuffmanNode* left = minHeap.pop();
        HuffmanNode* right = minHeap.pop();

        // 打印合并的节点信息
        std::cout << "Merging Nodes: (" << left->value << ", " << left->freq << ") and ("
                  << right->value << ", " << right->freq << ")\n";

        // 创建新的合并节点
        HuffmanNode* merged = new HuffmanNode(-1, left->freq + right->freq);
        merged->left = left;
        merged->right = right;

        // 打印合并后的新节点
        std::cout << "New Node: Value = " << merged->value << " (internal node), Frequency = " << merged->freq << "\n";

        // 将新的合并节点放回堆中
        minHeap.push(merged);

        // 打印当前堆的内容
        std::cout << "Current Heap:\n";
        std::vector<HuffmanNode*> tempHeap = minHeap.getHeap();  // 获取堆的副本
        for (HuffmanNode* node : tempHeap) {
            std::cout << "Node: Value = " << node->value << ", Frequency = " << node->freq << std::endl;
        }
        std::cout << "----------------------------------------\n";
    }

    // 返回哈夫曼树的根节点
    return minHeap.pop();
}

// 递归生成哈夫曼编码
void generateHuffmanCodes(HuffmanNode* root, const std::string& code,
                          std::unordered_map<int, std::string>& codes) {
    if (!root) return;

    // 叶子节点
    if (!root->left && !root->right) {
        codes[root->value] = code;
    }

    generateHuffmanCodes(root->left, code + "0", codes);
    generateHuffmanCodes(root->right, code + "1", codes);
}

// 保存哈夫曼编码表到文件
void saveHuffmanTableToFile(const std::unordered_map<int, std::string>& codes, const std::string& filename) {
    std::ofstream outFile(filename, std::ios::out);  // 打开文件以写入数据
    for (const auto& entry : codes) {
        outFile << entry.first << " " << entry.second << std::endl;  // 保存符号与其对应的编码
    }
    outFile.close();
}

// 从文件读取哈夫曼编码表
std::unordered_map<int, std::string> readHuffmanTableFromFile(const std::string& filename) {
    std::ifstream inFile(filename, std::ios::in);  // 打开文件以读取数据
    std::unordered_map<int, std::string> codes;
    int value;
    std::string code;
    while (inFile >> value >> code) {
        codes[value] = code;
    }
    inFile.close();
    return codes;
}

// 将哈夫曼编码序列存储到文件
void saveEncodedDataToFile(const std::vector<std::string>& encodedData, const std::string& filename) {
    std::ofstream outFile(filename, std::ios::out);  // 打开文件以写入数据
    for (const std::string& code : encodedData) {
        outFile << code << " ";  // 直接写入编码
    }
    outFile.close();
}

// 从文件读取编码数据
std::vector<std::string> readEncodedDataFromFile(const std::string& filename) {
    std::ifstream inFile(filename, std::ios::in);  // 打开文件以读取数据
    std::vector<std::string> encodedData;
    std::string code;
    while (inFile >> code) {
        encodedData.push_back(code);  // 按照空格分隔读取编码
    }
    inFile.close();
    return encodedData;
}

// 对数组进行哈夫曼编码压缩
std::vector<std::string> compressData(const std::vector<int>& data, std::unordered_map<int, std::string>& codes) {
    // 统计每个符号的频率
    std::unordered_map<int, int> freqMap;
    for (int value : data) {
        freqMap[value]++;
    }
    // 输出符号频率
    std::cout << "Symbol Frequencies:\n";
    for (const auto& entry : freqMap) {
        std::cout << "Value: " << entry.first << " Frequency: " << entry.second << std::endl;
    }


    // 构建哈夫曼树
    HuffmanNode* root = buildHuffmanTree(freqMap);

    // 生成哈夫曼编码
    generateHuffmanCodes(root, "", codes);

    // 将数据编码为哈夫曼编码
    std::vector<std::string> encodedData;
    for (int value : data) {
        encodedData.push_back(codes[value]);
    }

    // 将编码表和编码序列分别存储到文件
    saveHuffmanTableToFile(codes, "E:\\Codes\\Cpp\\Baselines\\pgm_index\\rsd_huffman\\huffman_table.txt");  // 存储哈夫曼编码表
    saveEncodedDataToFile(encodedData, "E:\\Codes\\Cpp\\Baselines\\pgm_index\\rsd_huffman\\encoded_data.txt");  // 存储编码序列

    return encodedData;
}

// 解码函数
std::vector<int> decodeData(const std::vector<std::string>& encodedData, const std::unordered_map<std::string, int>& reverseCodes) {
    std::vector<int> decodedData;
    std::string currentCode = "";

    for (const std::string& code : encodedData) {
        currentCode += code;  // 将每个编码片段添加到当前代码串

        if (reverseCodes.find(currentCode) != reverseCodes.end()) {
            decodedData.push_back(reverseCodes.at(currentCode));  // 找到匹配的编码
            currentCode = "";  // 重置当前代码以解码下一个符号
        }
    }

    return decodedData;
}

// 计算哈夫曼编码表和编码数据的总占用空间（以字节为单位）
void calculateSpaceUsage(const std::unordered_map<int, std::string>& codes, const std::vector<int>& data) {
    size_t totalTableSize = 0;  // 哈夫曼编码表占用的字节数
    size_t totalEncodedDataSize = 0;  // 编码数据占用的字节数

    // 计算哈夫曼编码表的大小
    for (const auto& entry : codes) {
        int symbolSize = sizeof(entry.first);  // 假设符号占用一个int的大小
        int codeSize = entry.second.size();  // 编码长度（单位是比特）
        totalTableSize += symbolSize + (codeSize / 8 + (codeSize % 8 ? 1 : 0));  // 每个编码项的总大小（单位字节）
    }

    // 计算编码数据的大小
    for (int value : data) {
        const std::string& code = codes.at(value);
        totalEncodedDataSize += (code.size() / 8 + (code.size() % 8 ? 1 : 0));  // 将每个符号的编码长度转为字节
    }

    // 输出结果
    std::cout << "Original Data Size: " << data.size() * sizeof(int) << " bytes" << std::endl;
    std::cout << "Huffman Table Size: " << totalTableSize << " bytes" << std::endl;
    std::cout << "Encoded Data Size: " << totalEncodedDataSize << " bytes" << std::endl;
    std::cout << "Total Space Usage: " << totalTableSize + totalEncodedDataSize << " bytes" << std::endl;
}

int main() {
    // 测试数据：假设我们有一个存储残差的数组，值都小于等于127
    std::vector<int> data = {5, 7, 5, 6, 7, 8, 5, 7, 7, 6, 8, 5, 10, 11, 12, 13, 14, 15, 127, 128, 129, 130, 5, 5, 5, 5, 5, 5, 5, 5, 5};

    // 保存编码表
    std::unordered_map<int, std::string> codes;

    // 对数据进行哈夫曼编码压缩
    std::vector<std::string> compressedData = compressData(data, codes);

    // 输出哈夫曼编码表
    // std::cout << "Huffman Encoding Table:\n";
    // for (const auto& entry : codes) {
    //     std::cout << "Value: " << entry.first << " Code: " << entry.second << std::endl;
    // }

    // // 从文件读取哈夫曼编码表并解码
    // std::unordered_map<int, std::string> readCodes = readHuffmanTableFromFile("E:\\Codes\\Cpp\\Baselines\\pgm_index\\rsd_huffman\\huffman_table.txt");
    //
    // // 反转编码表用于解码
    // std::unordered_map<std::string, int> reverseCodes;
    // for (const auto& entry : readCodes) {
    //     reverseCodes[entry.second] = entry.first;
    // }
    //
    // // 读取编码序列并解码
    // std::vector<std::string> encodedDataFromFile = readEncodedDataFromFile("E:\\Codes\\Cpp\\Baselines\\pgm_index\\rsd_huffman\\encoded_data.txt");
    // std::vector<int> decodedData = decodeData(encodedDataFromFile, reverseCodes);
    //
    // // 输出解压后的数据
    // std::cout << "\nDecoded Data: ";
    // for (int value : decodedData) {
    //     std::cout << value << " ";
    // }
    // std::cout << std::endl;

    calculateSpaceUsage(codes, data);

    return 0;
}
