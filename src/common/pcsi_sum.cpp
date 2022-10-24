
#include "pcsi_sum.h"

#include "osn.h"
#include "constants.h"

#include "poly/poly.h"
#include "ot/ot.h"

#include "ENCRYPTO_utils/connection.h"
#include "ENCRYPTO_utils/socket.h"

#include "HashingTables/cuckoo_hashing/cuckoo_hashing.h"
#include "HashingTables/simple_hashing/simple_hashing.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <unordered_set>
#include <vector>
#include <cmath>

#include <bits/stdc++.h>

namespace PCSI{

int debugx = 0;

std::unique_ptr<CSocket> create_socket(const std::string ip, uint16_t port, e_role role) {
    std::unique_ptr<CSocket> sock;
    if (role == SERVER) {
        sock = Listen(ip.c_str(), port);
    } else {
        sock = Connect(ip.c_str(), port);
    }
    return sock;
}

void interpolate_poly_with_dummy(std::vector<uint64_t>::iterator poly_offset, 
                 std::vector<uint64_t>::const_iterator x_offset, 
                 std::vector<std::vector<uint64_t>>::const_iterator y_offset,
                 std::size_t bin_num,
                 PCSIContext& ctx) {
    osuCrypto::PRNG prng(_mm_set_epi32(0, 0, 0, 2));
    std::vector<ZpLongEle> X(ctx.poly_size), Y(ctx.poly_size), co(ctx.poly_size);
    printf("X size = %d\n", X.size());
    for (auto i = 0ull, bin_index = 0ull; i < ctx.poly_size;) {
        if (bin_index < bin_num) {
            if ((*y_offset).size() > 0) {
                for (auto &mask : *y_offset) {
                    // X[i].ele = mask & _61_mask;
                    printf("i = %d\n", i);

                    X.at(i).ele  = mask & _61_mask;
                    Y.at(i).ele = X[i].ele ^ *(x_offset);
                    i++;
                    if (i > (ctx.poly_size-1)) break;
                }
            }
            x_offset++;
            y_offset++;
            bin_index++;
        } else {
            // X[i].ele = prng.get<uint64_t>();
            // Y[i].ele = prng.get<uint64_t>();
            // i++;
        }
    }

    // Poly::interpolate(co, X, Y);
    // auto coeff = co.begin();
    // for (auto i = 0ull; i < co.size(); i++, poly_offset++, coeff++) {
    //     *poly_offset = (*coeff).ele;
    // }

}

void interpolate_poly(std::vector<uint64_t>& polys, std::vector<uint64_t>& X, std::vector<std::vector<uint64_t>>& Y, PCSIContext& ctx) {
    std::size_t bins_num = Y.size();
    std::size_t Y_offset = 0;
    std::size_t bin_num_in_mega_bin = ceil_divide(bins_num, ctx.mega_bins_num);

    for (auto i = 0ull; i < ctx.mega_bins_num; i++) {
        auto poly = polys.begin() + ctx.poly_size * i;
        auto x = X.begin() + bin_num_in_mega_bin * i;
        auto y = Y.begin() + bin_num_in_mega_bin * i;

        if ((Y_offset + bin_num_in_mega_bin) > Y.size()) {
            auto overflow = (Y_offset + bin_num_in_mega_bin) % Y.size();
            bin_num_in_mega_bin -= overflow;
        }

        interpolate_poly_with_dummy(poly, x, y, bin_num_in_mega_bin, ctx);
        Y_offset += bin_num_in_mega_bin;
    }
}

void gen_corr_block(unsigned int n, int bins, int cur_level, int index, std::vector<uint64_t>& src, 
                    std::vector<std::array<osuCrypto::block, 2>>& ot_output, 
                    std::vector<osuCrypto::block>& corr_blocks) {
    // generate msg
    // msg[0] = m0 ^ w0 || m1 ^ w1
    // msg[1] = m0 ^ w1 || m1 ^ w0

    debugx += 1;
    int value_num = src.size();
   
    std::vector<uint64_t> bottom1;
    std::vector<uint64_t> top1;

    uint64_t m0, m1, w0, w1, tmp_m0[2], tmp_m1[2], corr_msg[2];
    osuCrypto::block tmp_block;


    if (value_num == 2) {
        if (n == 1) {
            m0 = src[0];
            m1 = src[1];
            tmp_block = ot_output[cur_level * (bins / 2) + index][0];
            memcpy(tmp_m0, &tmp_block, sizeof(tmp_m0));
            w0 = tmp_m0[0] ^ m0;
            w1 = tmp_m0[1] ^ m1;

            tmp_block = ot_output[cur_level * (bins / 2) + index][1];
            memcpy(tmp_m1, &tmp_block, sizeof(tmp_m1));

            corr_msg[0] = tmp_m1[0] ^ m0 ^ w1;
            corr_msg[1] = tmp_m1[1] ^ m1 ^ w0;

            corr_blocks[cur_level * (bins / 2) + index] = osuCrypto::toBlock(corr_msg[1], corr_msg[0]);

            tmp_m1[0] = m0 ^ w1;
            tmp_m1[1] = m1 ^ w0;
            ot_output[cur_level * (bins / 2) + index][1] = osuCrypto::toBlock(tmp_m1[1], tmp_m1[0]);
            src[0] = w0;
            src[1] = w1;
        } else {
            m0 = src[0];
            m1 = src[1];
            tmp_block = ot_output[(cur_level + 1) * (bins / 2) + index][0];
            memcpy(tmp_m0, &tmp_block, sizeof(tmp_m0));
            w0 = tmp_m0[0] ^ m0;
            w1 = tmp_m0[1] ^ m1;

            tmp_block = ot_output[(cur_level + 1) * (bins / 2) + index][1];
            memcpy(tmp_m1, &tmp_block, sizeof(tmp_m1));

            corr_msg[0] = tmp_m1[0] ^ m0 ^ w1;
            corr_msg[1] = tmp_m1[1] ^ m1 ^ w0;

            corr_blocks[(cur_level + 1) * (bins / 2) + index] = osuCrypto::toBlock(corr_msg[1], corr_msg[0]);

            tmp_m1[0] = m0 ^ w1;
            tmp_m1[1] = m1 ^ w0;
            ot_output[(cur_level + 1) * (bins / 2) + index][1] = osuCrypto::toBlock(tmp_m1[1], tmp_m1[0]);
            src[0] = w0;
            src[1] = w1;
        }
        return;
    }

    if (value_num == 3) {
        // process 1 and 2
        m0 = src[0];
        m1 = src[1];
        tmp_block = ot_output[cur_level * (bins / 2) + index][0];
        memcpy(tmp_m0, &tmp_block, sizeof(tmp_m0));
        w0 = tmp_m0[0] ^ m0;
        w1 = tmp_m0[1] ^ m1;

        tmp_block = ot_output[cur_level * (bins / 2) + index][1];
        memcpy(tmp_m1, &tmp_block, sizeof(tmp_m1));

        corr_msg[0] = tmp_m1[0] ^ m0 ^ w1;
        corr_msg[1] = tmp_m1[1] ^ m1 ^ w0;

        corr_blocks[cur_level * (bins / 2) + index] = osuCrypto::toBlock(corr_msg[1], corr_msg[0]);

        tmp_m1[0] = m0 ^ w1;
        tmp_m1[1] = m1 ^ w0;
        ot_output[cur_level * (bins / 2) + index][1] = osuCrypto::toBlock(tmp_m1[1], tmp_m1[0]);
        src[0] = w0;
        src[1] = w1;

        // process 2 and 3
        m0 = src[1];
        m1 = src[2];
        tmp_block = ot_output[(cur_level + 1) * (bins / 2) + index][0];
        memcpy(tmp_m0, &tmp_block, sizeof(tmp_m0));
        w0 = tmp_m0[0] ^ m0;
        w1 = tmp_m0[1] ^ m1;

        tmp_block = ot_output[(cur_level + 1) * (bins / 2) + index][1];
        memcpy(tmp_m1, &tmp_block, sizeof(tmp_m1));

        corr_msg[0] = tmp_m1[0] ^ m0 ^ w1;
        corr_msg[1] = tmp_m1[1] ^ m1 ^ w0;

        corr_blocks[(cur_level + 1) * (bins / 2) + index] = osuCrypto::toBlock(corr_msg[1], corr_msg[0]);
        
        tmp_m1[0] = m0 ^ w1;
        tmp_m1[1] = m1 ^ w0;
        ot_output[(cur_level + 1) * (bins / 2) + index][1] = osuCrypto::toBlock(tmp_m1[1], tmp_m1[0]);
        src[0] = w0;
        src[1] = w1;
        // process 1 and 3
        m0 = src[1];
        m1 = src[2];
        tmp_block = ot_output[(cur_level + 2) * (bins / 2) + index][0];
        memcpy(tmp_m0, &tmp_block, sizeof(tmp_m0));
        w0 = tmp_m0[0] ^ m0;
        w1 = tmp_m0[1] ^ m1;

        tmp_block = ot_output[(cur_level + 2) * (bins / 2) + index][1];
        memcpy(tmp_m1, &tmp_block, sizeof(tmp_m1));

        corr_msg[0] = tmp_m1[0] ^ m0 ^ w1;
        corr_msg[1] = tmp_m1[1] ^ m1 ^ w0;

        corr_blocks[(cur_level + 2) * (bins / 2) + index] = osuCrypto::toBlock(corr_msg[1], corr_msg[0]);
        
        tmp_m1[0] = m0 ^ w1;
        tmp_m1[1] = m1 ^ w0;
        ot_output[(cur_level + 2) * (bins / 2) + index][1] = osuCrypto::toBlock(tmp_m1[1], tmp_m1[0]);
        src[0] = w0;
        src[1] = w1;
        
        return;
    }
    
    int levels = 2 * n - 1;
    for (int i = 0; i < value_num - 1; i++) {
        m0 = src[i];
        m1 = src[i ^ 1];
        tmp_block = ot_output[cur_level * (bins / 2) + (i / 2)][0];
        memcpy(tmp_m0, &tmp_block, sizeof(tmp_m0));
        w0 = tmp_m0[0] ^ m0;
        w1 = tmp_m0[1] ^ m1;

        tmp_block = ot_output[cur_level * (bins / 2) + index][1];
        memcpy(tmp_m1, &tmp_block, sizeof(tmp_m1));

        corr_msg[0] = tmp_m1[0] ^ m0 ^ w1;
        corr_msg[1] = tmp_m1[1] ^ m1 ^ w0;

        corr_blocks[cur_level * (bins / 2) + index] = osuCrypto::toBlock(corr_msg[1], corr_msg[0]);
        
        tmp_m1[0] = m0 ^ w1;
        tmp_m1[1] = m1 ^ w0;
        ot_output[cur_level * (bins / 2) + index][1] = osuCrypto::toBlock(tmp_m1[1], tmp_m1[0]);
        src[i] = w0;
        src[i ^ 1] = w1;
        bottom1.push_back(src[i]);
        top1.push_back(src[i ^ 1]);
    }

    if (value_num & 1) {
        top1.push_back(src[value_num - 1]);
    }

    if (n > 0) gen_corr_block(n - 1, bins, cur_level + 1, index, bottom1, ot_output, corr_blocks);
    if (n > 0) gen_corr_block(n - 1, bins, cur_level + 1, index + value_num / 4, bottom1, ot_output, corr_blocks);

    for (int i = 0; i < value_num - 1; i+=2) {
        m0 = top1[i / 2];
        m1 = bottom1[i / 2];
        tmp_block = ot_output[(cur_level + levels - 1) * (bins / 2) + index + (i / 2)][0];
        memcpy(tmp_m0, &tmp_block, sizeof(tmp_m0));
        w0 = tmp_m0[0] ^ m0;
        w1 = tmp_m0[1] ^ m1;

        tmp_block = ot_output[(cur_level + levels - 1) * (bins / 2) + index + (i / 2)][1];
        memcpy(tmp_m1, &tmp_block, sizeof(tmp_m1));

        corr_msg[0] = tmp_m1[0] ^ m0 ^ w1;
        corr_msg[1] = tmp_m1[1] ^ m1 ^ w0;

        corr_blocks[(cur_level + levels - 1) * (bins / 2) + index + (i / 2)] = osuCrypto::toBlock(corr_msg[1], corr_msg[0]);
        
        tmp_m1[0] = m0 ^ w1;
        tmp_m1[1] = m1 ^ w0;
        ot_output[(cur_level + levels - 1) * (bins / 2) + index + (i / 2)][1] =osuCrypto::toBlock(tmp_m1[1], tmp_m1[0]);
        src[i] = w0;
        src[i ^ 1] = w1;
    }

    int middle = value_num >> 1;
    if (value_num & 1) {
        src[value_num - 1] = top1[middle];
    }
    // printf("leave gen_corr_block, value_num = %d, n = %d, debugx = %d , bottom1 = %d,\n", value_num, n, debugx, bottom1.size());
}

std::vector<osuCrypto::block> recv_osn(std::vector<int>& dest, int bins, PCSIContext& ctx) {
    // F_osn, alice as receiver, with input Pi
    std::vector<int> src(bins);
    for (int i = 0; i < bins; i++) {
        src[i] = i;
        dest[i] = i;
    }

    osuCrypto::PRNG prng(_mm_set_epi32(0, 0, 0, 0));

    // generate random permutation
    for (int i = bins - 1; i >= 0; i--) {
        int index = prng.get<uint64_t>() % (i + 1);
        dest[i] ^= dest[index];
        dest[index] ^= dest[i];
        dest[i] ^= dest[index];
    } 

    int n = int(ceil(log2(bins)));
    // generate transposition
    std::cout << "start route" << std::endl;
    route(n, 0, 0, src, dest);
    std::cout << "finish route" << std::endl;
    // the switch set is the choices in ot
    osuCrypto::BitVector choices = retrieve_switches(bins);
    std::cout << "finish get" << std::endl;

    // process ot
    osuCrypto::IOService ios;
    std::string name = "pcsi";
    std::vector<osuCrypto::block> recv_msg(choices.size());
    std::vector<osuCrypto::block> recv_corr(choices.size());

    rot_recv(choices, recv_msg, ctx);
    printf("finish rot_recv");

    osuCrypto::Session ep(ios, ctx.ip, ctx.port + 100, osuCrypto::SessionMode::Server, name);
    auto osn_recv_chl = ep.addChannel(name, name);

    osn_recv_chl.recv(recv_corr.data(), recv_corr.size()); // this statement error

    std::cout << "finish ot" << std::endl;
    // consider switch
    uint64_t tmp_msg[2], tmp_corr[2];
    for (int i = 0; i < recv_msg.size(); i++) {
        // choices means switches
        if (choices[i]) {
            tmp_msg[0] = (*((uint64_t *)&recv_msg)) ^ (*((uint64_t *)&recv_corr));
            tmp_msg[1] = (*((uint64_t *)&recv_msg) + 1) ^ (*((uint64_t *)&recv_corr) + 1);
            recv_msg[i] = osuCrypto::toBlock(tmp_msg[1], tmp_msg[0]);
        }
    }

    return recv_msg;
}

// SELECT SUM(A.VAL) FROM A UNION B
// SUM(A.VAL) + SUM(B.VAL) - SUM(A intersection B.VAL)

std::vector<std::vector<uint64_t>> send_osn(int bins, PCSIContext& ctx) {
    int n = int(ceil(log2(bins)));   // bins = 10, n = 4
    int levels = 2 * n - 1;   // 7
    int switch_num = levels * (bins / 2); // 7*5 = 35
    std::vector<uint64_t> masks(bins);
    std::vector<std::vector<uint64_t>> mat_masks(bins);


    osuCrypto::PRNG prng(_mm_set_epi32(0, 0, 0, 1));

    for (int i = 0; i < bins; i++) {
        uint64_t r = prng.get<uint64_t>();
        masks[i] = r;
        mat_masks[i].push_back(r);  // think;
    }

    std::vector<std::array<osuCrypto::block, 2>> ot_msg(switch_num);  // 电路中输入、输出、门的个数
    // ot_msg is valued in rot, it is outPut;
    rot_send(ot_msg, ctx);    // this is OT , ot_msg is output;
    printf("test send_osn,  switch_num = %d.\n", switch_num);

    std::vector<osuCrypto::block> corr_blocks(switch_num);
    
    std::cout << "n = " << n << std::endl;
    std::cout << "bins = " << bins << std::endl;
    gen_corr_block((unsigned int)n, bins, 0, 0, masks, ot_msg, corr_blocks);
    printf("leave gen_corr_block ********************************************** .\n");

    osuCrypto::IOService ios;
    std::string name = "pcsi";

    osuCrypto::Session ep(ios, ctx.ip, ctx.port + 100, osuCrypto::SessionMode::Client, name); 
    auto osn_send_chl = ep.addChannel(name, name);
    // osn_send_chl.asyncSend(corr_blocks);
    osn_send_chl.send(corr_blocks);

    for (int i = 0; i < bins; i++) {
        mat_masks[i].push_back(masks[i]);
    }
    return mat_masks;
}

std::vector<uint64_t> client_opprf(const std::vector<uint64_t>& eles, std::vector<uint64_t>& cuckoo_table, PCSIContext& ctx) {
    // 1. cuckoo hashing
    // initial cuckoo table
    ENCRYPTO::CuckooTable table(static_cast<size_t>(ctx.bins_num));
    table.SetNumOfHashFunctions(ctx.hash_num);
    table.Insert(eles);
    // map elements
    table.MapElements();
    // get vector
    cuckoo_table = table.AsRawVector();
    std::cout << "test, finish get vector ." << std::endl;

    // 2. oprf
    std::vector<uint64_t> masks = oprf_receiver(cuckoo_table, ctx);

    std::cout << "test, finish oprf receiver ." << std::endl;

    // 3. get polynominal and evaluate
    std::unique_ptr<CSocket> sock = create_socket(ctx.ip, ctx.port, static_cast<e_role>(ctx.role));
    const auto bin_num_in_mega_bin = ceil_divide(ctx.bins_num, ctx.mega_bins_num);
    std::vector<std::vector<ZpLongEle>> polys(ctx.mega_bins_num);
    std::vector<ZpLongEle> X(ctx.bins_num), Y(ctx.bins_num);
    for (auto &poly : polys) {
        poly.resize(ctx.poly_size);
    }

    for (auto i = 0ull; i < X.size(); i++) {
        X[i].ele = masks[i];
    }
    // 3.1 receive poly from server
    std::vector<uint8_t> poly_buffer(ctx.mega_bins_num * ctx.poly_bytelength);

    sock->Receive(poly_buffer.data(), ctx.mega_bins_num * ctx.poly_bytelength);
    sock->Close();
    std::cout << "test, get poly after socket recv ." << std::endl;

    // 3.2 evaluate
    for (auto i = 0; i < poly_buffer.size() && i < ctx.mega_bins_num; i++) {  // cyc:check
        for (auto j = 0; j < ctx.poly_size; j++) {
            // polys[i][j].ele = (reinterpret_cast<uint64_t *>(
            polys.at(i).at(j).ele = (reinterpret_cast<uint64_t *>(
          poly_buffer.data()))[i * ctx.poly_size + j];
        }
    }


    for (auto i = 0; i < X.size(); i++) {
        auto p = i / bin_num_in_mega_bin;
        Poly::eval(Y.at(i), polys.at(p), X.at(i));
    }

    // 3.3 xor
    std::vector<uint64_t> result;
    result.reserve(X.size());
    for (auto i = 0ull; i < X.size(); i++) {
        result.push_back(X[i].ele ^ Y[i].ele);
    }

    return result;
}

std::vector<uint64_t> server_opprf(const std::vector<uint64_t>& eles, PCSIContext& ctx) {

    ENCRYPTO::SimpleTable table(static_cast<size_t>(ctx.bins_num));
    table.SetNumOfHashFunctions(ctx.hash_num);
    table.Insert(eles);
    table.MapElements();

    auto simple_table = table.AsRaw2DVector();
    auto masks = oprf_sender(simple_table, ctx);

    std::vector<uint64_t> polys(ctx.mega_bins_num * ctx.poly_size, 0);
    std::vector<uint64_t> X(ctx.bins_num);

    std::random_device urandom("/dev/urandom");
    std::uniform_int_distribution<uint64_t> dist(1, (1ull << ctx.max_bitlen) - 1);
    std::generate(X.begin(), X.end(), [&]() { return dist(urandom); });
    auto tmp = X;
    std::sort(tmp.begin(), tmp.end());
    auto last = std::unique(tmp.begin(), tmp.end());
    tmp.erase(last, tmp.end());

    interpolate_poly(polys, X, masks, ctx);   // test is correctness.

    for (int xx = 0; xx < polys.size();xx++) {
        printf("%08x,", polys[xx]);
    }
    printf("\nploy vector print\n");

    std::unique_ptr<CSocket> sock = create_socket(ctx.ip, ctx.port, static_cast<e_role>(ctx.role));
    sock->Send((uint8_t *)polys.data(), ctx.mega_bins_num * ctx.poly_bytelength); // 已测试，不是sock的问题
    sock->Close();

    std::cout << "socket has setup" << std::endl;

    return X;
}

uint64_t exec(const std::vector<uint64_t>& inputs, PCSIContext& ctx, const std::vector<uint64_t>& data) {
    // network
    printf("enter exec.\n");
    std::unique_ptr<CSocket> sock = create_socket(ctx.ip, ctx.port, static_cast<e_role>(ctx.role));
    sock->Close();

    std::vector<uint64_t> input_bak(inputs), sets;
    int bins = ctx.bins_num;   // coocku hash bins num: hash值的长度？ = n_eles * 1.27；与元素个数线性关系
    std::cout << "ctx.bins_num = " << bins << std::endl;
    uint64_t output = 0;

    if (ctx.role == CLIENT) {
        // client give Data
        // 1. osn preprocessing 
        // 1.1 initial paramters
        std::vector<int> dest(bins);
        int n = int(ceil(log2(bins)));
        int levels = 2 * n - 1;
        initialize(bins, levels);
	    std::cout << "[client]finish 1.1" << std::endl;
        // 1.2 processing offline osn
        std::vector<osuCrypto::block> ot_output;
        ot_output = recv_osn(dest, bins, ctx);
        // dest : x_Pi(w) ^ Mask, is output;
	    std::cout << "[client]finish 1" << std::endl;
        // 2. pcsi preprocessing 
        std::vector<uint64_t> cuckoo_table;
        sets = client_opprf(input_bak, cuckoo_table, ctx);
	    std::cout << "[client]finish 2" << std::endl;
        // // 3. online osn
        // std::vector<uint64_t> shuffled_sets(bins);
        // for (int i = 0; i < bins; i++) {
        //     shuffled_sets[i] = sets[dest[i]];
        // }

        // std::string name = "pcsi";
        // osuCrypto::IOService ios;
        // osuCrypto::Session ep(ios, ctx.ip, ctx.port + 20, osuCrypto::SessionMode::Server, name);
        // auto recv_chl = ep.addChannel(name, name);
        // std::vector<uint64_t> input_vec(bins);
        // recv_chl.recv(input_vec.data(), input_vec.size());

        // std::vector<std::vector<osuCrypto::block>> mat_ot_output(levels, std::vector<osuCrypto::block>(bins));
        // int index = 0;
        // for (int i = 0; i < levels; i++) {
        //     for (int j = 0; j < bins / 2; j++) {
        //         mat_ot_output[i][j] = ot_output[index++];
        //     }
        // }
        // masked_eval(n, 0, 0, input_vec, mat_ot_output);
	    // std::cout << "[client]finish 3" << std::endl;
        // // 4. send the result to process equality test
        // std::string name_ = "pcsi";
        // osuCrypto::Session ep_(ios, ctx.ip, ctx.port + 30, osuCrypto::SessionMode::Server, name_);
        // auto send_chl_eq = ep_.addChannel(name_, name_);
        // std::vector<uint64_t> sets_eq;
        // for (int i = 0; i < bins; i++) {
        //     sets_eq.push_back(shuffled_sets[i] ^ input_vec[i]);
        // }
        // send_chl_eq.asyncSend(sets_eq);
	    // std::cout << "[client]finish 4" << std::endl;
        // // 5. do sum
        // std::vector<std::pair<uint64_t, uint64_t>> candidate;
        // candidate.reserve(input_bak.size());
        // std::transform(input_bak.begin(), input_bak.end(), data.begin(), std::back_inserter(candidate), [](uint64_t x, uint64_t y){ return std::make_pair(x, y); });
        // std::sort(candidate.begin(), candidate.end());

        // std::vector<uint64_t> candidate_index(candidate.size());
        // for (int i = 0; i < candidate.size(); i++) {
        //     candidate_index[i] = candidate[i].first;
        // }

        // // assign 
        // for (int i = 0; i < cuckoo_table.size(); i++) {
        //     std::vector<uint64_t>::iterator it = std::lower_bound(candidate_index.begin(), candidate_index.end(), cuckoo_table[i]);

        //     if (it == std::end(candidate_index)) {
        //         cuckoo_table[i] = 0;
        //     } else {
        //         cuckoo_table[i] = candidate[it-candidate_index.begin()].second;
        //     }
        // }

        // std::vector<uint64_t> shuffled_table(cuckoo_table.size());
        // for (int i = 0; i < cuckoo_table.size(); i++) {
        //     shuffled_table[i] = cuckoo_table[dest[i]];
        // }

        // osuCrypto::PRNG prng(_mm_set_epi32(0, 0, 0, 0));
        // std::vector<std::vector<osuCrypto::block>> msg(cuckoo_table.size());

        // for (int i = 0; i < cuckoo_table.size(); i++) {
        //     int r = prng.get<uint64_t>();
        //     msg[i].push_back(osuCrypto::toBlock(0, r));
        //     msg[i].push_back(osuCrypto::toBlock(0, shuffled_table[i] - r));
        //     output += r;
        // }
        // ot_send(msg, ctx);

    } else {
        // offline osn
        // server is receiver, have func, not data;
        std::cout << "server inter exec. " << std::endl;
        std::vector<std::vector<uint64_t>> pre_masks;
        pre_masks = send_osn(bins, ctx);
        printf("server send_osn, pre_masks size = %d .\n", pre_masks.size());

        // pcsi preprocessing
        sets = server_opprf(input_bak, ctx);

        printf("finish opprf.\n");

        // // online osn
        // osuCrypto::IOService ios;
        // std::string name = "pcsi";
        // osuCrypto::Session ep(ios, ctx.ip, ctx.port + 20, osuCrypto::SessionMode::Client, name);
        // auto send_chl = ep.addChannel(name, name);
        // std::vector<uint64_t> output_masks, benes_input;

        // for (int i = 0; i < sets.size(); i++) {
        //     pre_masks[i][0] ^= sets[i];
        // }
        // for (int i = 0; i < bins; i++) {
        //     benes_input.push_back(pre_masks[i][0]);
        // }
        // send_chl.asyncSend(benes_input);
        // for (int i = 0; i < ctx.bins_num; i++) {
        //     output_masks.push_back(pre_masks[i][1]);
        // }
        // printf("finish masks .\n");
        // // equality test
        // std::string name_ = "pcsi";
        // osuCrypto::Session ep_(ios, ctx.ip, ctx.port + 30, osuCrypto::SessionMode::Client, name_);
        // auto recv_chl_eq = ep_.addChannel(name_, name_);
        // std::vector<uint64_t> sets_eq(bins);
        // recv_chl_eq.recv(sets_eq.data(), sets_eq.size());
        // osuCrypto::BitVector char_vec(sets_eq.size());
        // for (int i=0; i < sets_eq.size();++i) {
        //     if (sets_eq[i] == output_masks[i]) {
        //         char_vec[i] = 1;
        //     }
        // }

        // // do sum
        // printf("begin do sum.\n");
        // std::vector<osuCrypto::block> recv_msg(char_vec.size());
        // ot_recv(char_vec, recv_msg, ctx);
        // uint64_t ot_msg[2];
        // for (int i=0; i < recv_msg.size(); ++i) {
        //     memcpy(ot_msg, &recv_msg[i], sizeof(ot_msg)); 
        //     output += ot_msg[0];
        // }
    }
    return output;
}

}

