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

#include "utils.hpp"

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_bitcompressed.hpp>
#include <sdsl/rmq_support.hpp>

#include <string_view>
#include <map>

template<typename CharT>
struct LZEndPhrase {
    size_t length;       ///< The length of the phrase (which includes the trailing character).
    size_t source;       ///< The phrase identifier where the source of this phrase ends.
    CharT trailing_char; ///< The trailing character of the phrase.
    LZEndPhrase(size_t length, size_t source, CharT trailing_char)
        : length(length), source(source), trailing_char(trailing_char) {}
};

struct LZEndParsingStrategyKK;
struct LZEndParsingStrategyKN;

/**
 * A class that constructs and stores in memory the LZ-End parsing of a given text.
 *
 * @tparam CharT character type
 */
template<typename CharT, typename ParsingStrategy = LZEndParsingStrategyKN>
class LZEnd {
    std::vector<LZEndPhrase<CharT>> phrases;
    size_t text_length;

public:

    LZEnd() = default;

    template<class Text>
    explicit LZEnd(Text &&s) : text_length(s.size()) {
        static_assert(std::is_same_v<CharT, typename std::remove_reference_t<Text>::value_type>,
                      "CharT must be the same as the value type of the text.");
        ParsingStrategy::parse(s, phrases);
    }

    /** Returns the decoded text. */
    sdsl::int_vector<sizeof(CharT) * 8> decode() const {
        sdsl::int_vector<sizeof(CharT) * 8> result(text_length + 1);
        sdsl::int_vector<> phrase_end(text_length, 0, sdsl::bits::hi(text_length) + 1);

        size_t i = 0;
        for (auto j = 0; j < phrases.size(); ++j) {
            auto &p = phrases[j];
            if (p.length > 1)
                std::copy_n(result.begin() + phrase_end[p.source] - p.length + 2, p.length - 1, result.begin() + i);
            i += p.length;
            result[i - 1] = p.trailing_char;
            phrase_end[j] = i - 1;
        }

        result.resize(text_length);
        return result;
    }

    /** Pretty-prints the LZ-End parsing of the text, using | as the phrase delimiter. */
    friend std::ostream &operator<<(std::ostream &os, const LZEnd<CharT> &parsing) {
        auto s = parsing.decode();

        size_t i = 0;
        for (auto &p: parsing.phrases) {
            std::copy_n(s.begin() + i, p.length, std::ostream_iterator<CharT>(os));
            i += p.length;
            if (i < parsing.data_size())
                os << "|";
        }

        return os;
    }

    auto data_size() const { return text_length; }
    auto size() const { return phrases.size(); }
    auto begin() const { return phrases.begin(); }
    auto end() const { return phrases.end(); }
    auto &operator[](size_t i) const { return phrases[i]; }
};

template<class Text>
auto lzend_parse(Text &&s) {
    return LZEnd<typename std::remove_reference_t<Text>::value_type>(std::forward<Text>(s));
}

template<class CharT>
auto lzend_parse(const CharT *s) { return LZEnd<CharT>(std::basic_string_view<CharT>(s)); }

/**
 * Implements the LZ-End parsing algorithm described in:
 *   Sebastian Kreft and Gonzalo Navarro. On compressing and indexing repetitive sequences. Theor. Comput. Sci. (2013).
 */
struct LZEndParsingStrategyKN {
    template<typename Text, typename CharT>
    static void parse(Text &&s, std::vector<LZEndPhrase<CharT>> &phrases) {
        using rk_t = sdsl::rank_support_v<1>;
        using sl1_t = sdsl::select_support_scan<1>;
        using sl0_t = sdsl::select_support_scan<0>;
        using csa_t = std::conditional_t<
            sizeof(CharT) == 1,
            sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector, rk_t, sl1_t, sl0_t>, 1, 1>,
            sdsl::csa_wt_int<sdsl::wt_huff_int<sdsl::bit_vector, rk_t, sl1_t, sl0_t>, 1, 1>>;
        csa_t sa;
        sdsl::rmq_support_sparse_table<decltype(sa), false> rmq;
        construct_im_reversed(s, sa, rmq);

        auto text_length = s.size();
        std::map<size_t, size_t> position_to_phrase; // map ending positions in the text to their phrase id
        position_to_phrase.emplace(text_length + 1, 0);
        size_t i = 0;
        while (i < text_length) {
            typename decltype(sa)::size_type sp = 0, ep = text_length;
            size_t ip = i;
            size_t j = i;
            size_t q = 0;

            while (ip < text_length) {
                auto matches = sdsl::backward_search(sa, sp, ep, s[ip], sp, ep);
                if (!matches || sa[rmq(sp, ep)] <= text_length - i - 1)
                    break;
                ++ip;
                auto[fpos, qp] = *position_to_phrase.lower_bound(sp);
                if (fpos <= ep) {
                    j = ip;
                    q = qp;
                }
            }

            if (j == text_length) {
                phrases.emplace_back(j - i + 1, q, CharT());
                break;
            }

            position_to_phrase.emplace(sa.isa[text_length - j - 1], phrases.size());
            phrases.emplace_back(j - i + 1, q, s[j]);
            i = j + 1;
        }
    }

    template<class Text, class SA, class RMQ>
    static void construct_im_reversed(Text &&s, SA &sa, RMQ &rmq) {
        using namespace sdsl;
        auto tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
        store_to_file_reversed(s, tmp_file);
        cache_config config(false, "@");
        construct(sa, tmp_file, config, 0);
        config.delete_files = true;
        config.delete_data = true;
        rmq = {&sa};
        ram_fs::remove(tmp_file);
    }
};

/**
 * Implements a data structure to be used in a LZ-End parsing of a given text.
 *
 * The data structure consists of a number of classical arrays built on the reversal of the given text, namely, the
 * suffix array, the inverse suffix array, the LCP array, and the RMQ data structure.
 *
 * Reference: Dominik Kempa and Dmitry Kosolobov. LZ-End Parsing in Linear Time. ESA 2017.
 *
 * @tparam CharT character type
 */
template<typename CharT>
class LCPStructure {
    sdsl::csa_bitcompressed<sdsl::null_alphabet<sdsl::int_alphabet_tag>> m_sa;
    sdsl::lcp_bitcompressed<> m_lcp;
    sdsl::rmq_support_sparse_table<decltype(m_lcp), true> m_rmq;
//     sdsl::rmq_succinct_sada<true> m_rmq;

public:

    /** Returns the ending position of the ith lexicographically smallest reverse prefix in the input text T[0,n-1],
     * i.e. sa() is such that reverse(T[0,sa(0)]) < reverse(T[0,sa(1)]) < ... < reverse(T[n-1,sa(n-1)]) */
    size_t sa(size_t i) const { return i == 0 ? m_sa[i] : m_sa.size() - 2 - m_sa[i]; }

    /** Returns the position i such that j == sa(i).  */
    size_t isa(size_t j) const { return j == m_sa.size() - 1 ? m_sa.isa[j] : m_sa.isa[m_sa.size() - 2 - j]; }

    /** Returns the length of the longest common suffix of T[0,i] and T[0,j]. */
    size_t longest_common_suffix(size_t i, size_t j) const {
        assert(i <= j);
        if (isa(i) > isa(j))
            return m_lcp[m_rmq(isa(j) + 1, isa(i))];
        return m_lcp[m_rmq(isa(i) + 1, isa(j))];
    }

    /** Constructs the LCP structure on a given object. The object must be a string or an int_vector. */
    template<class Text>
    explicit LCPStructure(Text &&text) {
        construct_im_reversed(text);

//        if constexpr (std::is_convertible_v<Text, std::string_view>) { // Print rotations
//            for (int i = 0; i <= text.size(); ++i) {
//                std::string a(text.substr(0, sa(i) + 1));
//                std::cout << i << ": " << a << " ";
//                if (i == 0) {
//                    std::cout << "$" << std::endl;
//                    continue;
//                }
//                std::reverse(a.begin(), a.end());
//                std::cout << a << std::endl;
//            }
//            for (int i = 0; i < text.size(); ++i) // Print all computations of longest_common_suffix
//                for (int j = i + 1; j < text.size(); ++j)
//                    std::cout << "i=" << i << " j=" << j
//                              << " isa(i)=" << isa(i) << " isa(j)=" << isa(j)
//                              << " " << text.substr(0, i + 1) << " " << text.substr(0, j + 1)
//                              << " LCSlen=" << longest_common_suffix(i, j) << std::endl;
//        }
    }

    /** Returns the size of the string, including the terminator. */
    size_t size() const { return m_sa.size(); }

private:

    template<class Text>
    void construct_im_reversed(Text &&text) {
        using namespace sdsl;
        auto tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
        store_to_file_reversed(text, tmp_file);
        cache_config config(false, "@");
        construct(m_sa, tmp_file, config, 0);
        config.delete_files = true;
        config.delete_data = true;
        construct(m_lcp, tmp_file, config, 0);
        m_rmq = {&m_lcp};
        ram_fs::remove(tmp_file);
    }
};

/**
 * Implements the LZ-End parsing algorithm described in:
 *   Dominik Kempa and Dmitry Kosolobov. LZ-End Parsing in Linear Time. ESA 2017. Section 2.
 */
struct LZEndParsingStrategyKK {
    template<typename Text, typename CharT>
    static void parse(Text &&s, std::vector<LZEndPhrase<CharT>> &phrases) {
        auto text_length = s.size();
        if (text_length < 3) {
            if (text_length >= 1) phrases.emplace_back(1, 0, s[0]);
            if (text_length >= 2) phrases.emplace_back(1, 0, s[1]);
            return;
        }

        LCPStructure<CharT> lcps(s);
        std::map<size_t, size_t> marked_prefixes; // map prefixes that end at phrase boundaries to their phrase id

        // Function to compute the predecessor and successor of a key in the set: marked_prefixes \ {ignore_key}
        auto nil = marked_prefixes.end();
        auto pred_succ = [&](size_t k, size_t ignore_key) {
            auto it = marked_prefixes.lower_bound(k);
            auto predecessor = it == marked_prefixes.begin() ? nil : std::prev(it);
            if (it != marked_prefixes.end() && it->first == k)
                ++it;
            auto successor = it == marked_prefixes.end() ? nil : it;
            if (predecessor != nil && predecessor->first == ignore_key)
                predecessor = predecessor == marked_prefixes.begin() ? nil : std::prev(predecessor);
            if (successor != nil && successor->first == ignore_key)
                successor = successor == marked_prefixes.end() ? nil : std::next(successor);
            return std::make_tuple(predecessor, successor);
        };

        // Function to find a candidate earlier occurrence of the last phrase in s[0,k] ending at a phrase boundary. It
        // returns the length and the id of the phrase where the candidate occurrence ends.
        auto earlier_occurrence = [&](size_t k, size_t ignore_key = std::numeric_limits<size_t>::max()) {
            auto isak = lcps.isa(k);
            auto[j_pred, j_succ] = pred_succ(isak, ignore_key);
            auto lcs1 = j_pred != nil ? lcps.longest_common_suffix(lcps.sa(j_pred->first), k) : 0;
            auto lcs2 = j_succ != nil ? lcps.longest_common_suffix(lcps.sa(j_succ->first), k) : 0;
            if (lcs1 > lcs2)
                return std::make_pair(lcs1, j_pred->second);
            return std::make_pair(lcs2, j_succ->second);
        };

        // Compute first two phrases
        size_t start;
        if (s[0] != s[1]) {
            phrases.emplace_back(1, 0, s[0]);
            phrases.emplace_back(1, 0, s[1]);
            marked_prefixes.emplace(lcps.isa(0), 0);
            marked_prefixes.emplace(lcps.isa(1), 1);
            start = 1;
        } else {
            phrases.emplace_back(1, 0, s[0]);
            phrases.emplace_back(2, 0, s[2]);
            marked_prefixes.emplace(lcps.isa(0), 0);
            marked_prefixes.emplace(lcps.isa(2), 1);
            start = 2;
        }

        // Compute the rest of the phrases
        for (size_t k = start; k < s.size(); ++k) {
            // s[0,k] has already been parsed
            auto new_char = k == s.size() - 1 ? CharT() : s[k + 1];

            // Check whether the last two phrases have an earlier occurrence in the text
            auto &last_phrase = phrases.back();
            auto &penultimate_phrase = phrases[phrases.size() - 2];
            auto ignored_mark = lcps.isa(k - last_phrase.length);
            auto[length, source_id] = earlier_occurrence(k, ignored_mark);
            if (length >= last_phrase.length + penultimate_phrase.length) {
                // Unite the new character with the last two phrases
                marked_prefixes.erase(ignored_mark);
                marked_prefixes.erase(lcps.isa(k));
                marked_prefixes.emplace(lcps.isa(k + 1), phrases.size() - 2);
                penultimate_phrase.source = source_id;
                penultimate_phrase.length += last_phrase.length + 1;
                penultimate_phrase.trailing_char = new_char;
                phrases.pop_back();
                continue;
            }

            // Check whether the last phrase has an earlier occurrence in the text
            std::tie(length, source_id) = earlier_occurrence(k);
            if (length >= last_phrase.length) {
                // Unite the new character with the last phrase
                marked_prefixes.erase(lcps.isa(k));
                marked_prefixes.emplace(lcps.isa(k + 1), phrases.size() - 1);
                ++last_phrase.length;
                last_phrase.source = source_id;
                last_phrase.trailing_char = new_char;
                continue;
            }

            // Add a new phrase with a single character
            phrases.emplace_back(1, 0, new_char);
            marked_prefixes.emplace(lcps.isa(k + 1), phrases.size() - 1);
        }
        if (phrases.back().length == 1 && phrases.back().trailing_char == CharT())
            phrases.pop_back();
    }
};