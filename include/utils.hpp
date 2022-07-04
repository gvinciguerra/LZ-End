// This file is part of LZ-End <https://github.com/gvinciguerra/LZ-End>.
// Copyright (c) 2021 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <sdsl/int_vector.hpp>

#include <type_traits>

namespace sdsl {

/** An alphabet strategy that provide no mapping. */
template<typename t_alphabet_tag>
struct null_alphabet {
    struct null_container;
    static constexpr bool is_byte_alphabet = std::is_same_v<t_alphabet_tag, sdsl::byte_alphabet_tag>;

    typedef int_vector<>::size_type size_type;
    typedef null_container char2comp_type;
    typedef null_container comp2char_type;
    typedef null_container C_type;
    typedef std::conditional_t<is_byte_alphabet, uint16_t, uint64_t> sigma_type;
    typedef std::conditional_t<is_byte_alphabet, uint8_t, uint64_t> char_type;
    typedef std::conditional_t<is_byte_alphabet, uint8_t, uint64_t> comp_char_type;
    typedef std::vector<char_type> string_type;
    typedef t_alphabet_tag alphabet_category;
    enum {
        int_width = t_alphabet_tag::WIDTH
    };

    struct null_container {
        char_type operator[](char_type const) const { throw std::logic_error("not implemented"); }
    };

    char2comp_type char2comp;
    comp2char_type comp2char;
    C_type C;
    sigma_type sigma;

    null_alphabet() = default;

    template<uint8_t t_width>
    null_alphabet(int_vector_buffer<t_width> &, int_vector_size_type) {}

    size_type serialize(std::ostream &, structure_tree_node *v, std::string const &name = "") const { return 0; }

    void load(std::istream &) {}
};

bool store_to_file_reversed(std::string_view v, const std::string &file) {
    using namespace sdsl;
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        std::cerr << "ERROR: store_to_file_reversed(const char *v, const std::string&)" << std::endl;
        return false;
    }

    sdsl::int_vector<8>::write_header(v.size() * 8, 8, out);
    std::copy(v.rbegin(), v.rend(), std::ostream_iterator<char>(out));
    out.close();
    return true;
}

template<uint8_t t_width>
bool store_to_file_reversed(const sdsl::int_vector<t_width> &v, const std::string &file) {
    using namespace sdsl;
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        std::cerr << "ERROR: store_to_file_reversed: Could not open file `" << file << "`" << std::endl;
        return false;
    }

    sdsl::int_vector<t_width> v_reversed = v;
    std::reverse(v_reversed.begin(), v_reversed.end());
    v_reversed.serialize(out, nullptr, "");
    out.close();
    return true;
}

}