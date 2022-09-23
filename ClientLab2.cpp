/*
 * File:   Main.cpp
 * Author: maslov_a
 *
 * Created on 17 сентября 2022 г., 13:56
 */

// g++ ClientLab2.cpp -o ClientLab2 -lpthread -lboost_program_options -fopenmp

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

#include <boost/asio.hpp>                   // TCP
#include <boost/numeric/ublas/matrix.hpp>   // Matrix
#include <boost/array.hpp>                  // Matrix
#include <boost/program_options.hpp>        // program args

// OpenMP
#include <omp.h>


using namespace std;
namespace basio  = boost::asio;
namespace bip    = boost::asio::ip;
namespace buplas = boost::numeric::ublas;
namespace po     = boost::program_options;

bool PROGRESS_BAR = false;

string gen_answer(vector<vector<double>> result_vector,
        chrono::high_resolution_clock::duration &duration) {
    string answer_str;

    // start measure
    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

    cout << "[Client] Generate answer ..." << endl;
    for (int i = 0; i < result_vector.size(); i++) {
        // header of matrix
        answer_str+=("Matrix[" + to_string(i) + "]\n");
        for (int j = 0; j < result_vector[i].size(); j++) {
            answer_str+=(to_string(result_vector[i][j]) + '\n');
        }
    }

    // stop measure
    chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    return answer_str;
}

vector<vector<double>> get_max(vector<buplas::matrix<double> > &matrix_vector,
        chrono::high_resolution_clock::duration &duration){

    cout << "[Client] Calc max numbers ..." << endl;

    // start measure
    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
    int vector_size = matrix_vector.size();
    vector<vector<double>> answer(vector_size);
    int matrix_size = matrix_vector[0].size1();

    #pragma omp parallel for
    for (int i = 0; i < vector_size; i++) {
        vector<double> max_lines(matrix_size);
        #pragma omp parallel for num_threads(matrix_size)
        for (int h = 0; h < matrix_size; h++) {
            vector<double> line(matrix_size);
            #pragma omp parallel for num_threads(matrix_size)
            for (int l = 0; l < matrix_size; l++) {
                line[l] = matrix_vector[i](h, l);
            }

            double max = *max_element(line.begin(), line.end());
            double result = max * max;
            max_lines[h] = result;
        }

        answer[i] = max_lines;
        // progress bar
        if (PROGRESS_BAR)
            printf("%d matrices of %d\r", i, vector_size);
    }

    // stop measure
    chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << endl;

    return answer;
}


// make vector of matrix from string
vector<buplas::matrix<double> >parse_input (string input,
        chrono::high_resolution_clock::duration &duration) {
    vector< buplas::matrix<double> > matrix_vector;

    // find first line (size of matrix)
    int size_pos = input.find('\n');
    int size     = stoi(input.substr(0, size_pos));
    int pos      = size_pos + 1;
    cout << "[Client] Matrix size = " << size <<  endl;

    int h, l;                                   // current coords of matrix
    int num_pos;                                // number position
    string num_str;                             // found number string
    double num;                                 // found number

    // start measure
    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
    int input_size = input.size();

    cout << "[Client] Parse input ... "  <<  endl;
    while (pos < input_size) {
        // create matrix
        buplas::matrix<double> matrix(size,size);

        // TODO: back old algorythm
        // parse height
        h = 0;
        while ( h < size ) {
            // parse length
            l = 0;
            while ( l < size ) {
                num_pos = input.find_first_of(" \n\0", pos);
                if (num_pos <= 0)
                    num_pos = input.size();
                num_str = input.substr(pos, num_pos - pos);
                num = stod(num_str);
                pos = num_pos + 1;

                matrix.insert_element(h, l , num);
                l++;
            }
            h++;
        }

        // add matrix
        matrix_vector.push_back(matrix);

        // progress bar
        if (PROGRESS_BAR)
            printf("%d bytes of %d\r", pos, input_size);
    }

    // stop measure
    chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << endl;

    return matrix_vector;
}


void print_matrix_vector(vector<buplas::matrix<double> > &matrix_vector) {
    int i;
    cout << endl;
    for (i = 0; i < matrix_vector.size(); i++){
        cout << "Matrix[" << i << "]" << endl;
        for (int h = 0; h < matrix_vector[i].size1(); h++) {
            for (int l = 0; l < matrix_vector[i].size2(); l++) {
                cout << matrix_vector[i](h, l) << " ";
            }
            cout << endl;
        }
    }
    cout << endl;
}

// make string from input tcp buffer
string make_string(basio::streambuf& streambuf) {
    return {buffers_begin(streambuf.data()), buffers_end(streambuf.data())};
}

// process args
int get_args(po::options_description &desc, int &P, string &H,
                                            int &ac, char **av) {
    desc.add_options()
        ("help,h", "print help")
        ("port,p", po::value<int>(),    "set port [default: 55555]")
        ("host,H", po::value<string>(), "set host ip [default: 127.0.0.1]")
        ("progress,P", po::value<bool>(), "progress bar 1/0 [default: 1]");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        exit(1);
    }

    if (vm.count("port"))
        P = vm["port"].as<int>();

    if (vm.count("host"))
        H = vm["host"].as<string>();

    if (vm.count("progress"))
        PROGRESS_BAR = vm["progress"].as<bool>();

    return 0;
}


int main(int argc, char* argv[]) {
    // GLOBAL VARS
    int PORT    = 55555;
    string HOST = "127.0.0.1";

    // times
    chrono::high_resolution_clock::duration parse_time, gen_answer_time, calcs_time;
    long double parse_time_s, gen_answer_time_s, calcs_time_s;

    // process args
    po::options_description desc("Options");
    get_args(desc, PORT, HOST, argc, argv);

    // Connect backend
    basio::io_context io_context;           // core IO object
    boost::system::error_code error;        // contain error
    bip::tcp::endpoint ep(bip::address::from_string(HOST), PORT);
    bip::tcp::socket socket(io_context);

    // try to connect to socket
    try {
        socket.connect(ep);
    } catch (const boost::system::system_error& ex) {
        cerr << "[Client error] Start server at first!" << endl;
        exit(1);
    }

    // Size of TCP message
    int HEADER_SIZE = 8192;
    int ACK_SIZE = 3;
    vector<char> buf(HEADER_SIZE);

    // RECEIVE
    // -------------------------------------------------------------------------
    cout << "[Client] Recieve message size ..." <<  endl;
    socket.read_some(basio::buffer(buf, HEADER_SIZE), error);
    if (error.value() != boost::system::errc::success)
        throw boost::system::system_error(error);
    int msg_size = stoi(buf.data());
    cout << "[Client] Message size: " << msg_size << endl;

    cout << "[Client] Send ack ..." << endl;
    basio::write(socket, basio::buffer("ACK", ACK_SIZE), error);

    buf.resize(msg_size);
    cout << "[Client] Recieve message ..." <<  endl;
    basio::read(socket, basio::buffer(buf, msg_size), error);
    if (error.value() != boost::system::errc::success)
        throw boost::system::system_error(error);
    string input(buf.begin(), buf.end());
    // -------------------------------------------------------------------------

    // CALCS
    // -------------------------------------------------------------------------
    // parse matrices
    cout << "[Client] Stream size: " << input.size() << " bytes" << endl;
    vector<buplas::matrix<double> > matrix_vector = parse_input(input, parse_time);
    cout << "[Client] Matrices number = " << matrix_vector.size() << endl;
    // time
    parse_time_s = parse_time.count()*1e-9;
    cout << "Parse time: " << parse_time_s << " s" << endl;

    // main buisness calcs
    vector<vector<double>> result_vector = get_max(matrix_vector, calcs_time);
    calcs_time_s = calcs_time.count()*1e-9;
    cout << "Calcs time: " << calcs_time_s << " s" << endl;

    // gen string from answer
    string ANSWER = gen_answer(result_vector, gen_answer_time);
    gen_answer_time_s = gen_answer_time.count()*1e-9;
    cout << "Gen answer time: " << gen_answer_time_s << " s"  << endl;

    // add times to answer
    string parse_time_str      = "Parse time, s:      ";
    string calcs_time_str      = "Calcs time, s:      ";
    string gen_answer_time_str = "Gen answer time, s: ";
    parse_time_str      += (to_string(parse_time_s) + "\n");
    calcs_time_str      += (to_string(calcs_time_s) + "\n");
    gen_answer_time_str += (to_string(gen_answer_time_s) + "\n");
    ANSWER += (parse_time_str + calcs_time_str + gen_answer_time_str);
    // -------------------------------------------------------------------------

    // SEND
    // -------------------------------------------------------------------------
    msg_size = ANSWER.size();
    cout << "[Client] Send message size ..." << endl;
    string msg_size_str = to_string(msg_size);
    basio::write(socket, basio::buffer(msg_size_str, HEADER_SIZE), error);

    buf.resize(HEADER_SIZE);
    cout << "[Client] Receive ack ..." <<  endl;
    basio::read(socket, basio::buffer(buf, ACK_SIZE), error);
    if (error.value() != boost::system::errc::success)
        throw boost::system::system_error(error);

    cout << "[Client] Send message ..." << endl;
    basio::write(socket, basio::buffer(ANSWER, msg_size), error);
    // -------------------------------------------------------------------------

    cout << "[Client] Finish" << endl;
    return 0;
}

