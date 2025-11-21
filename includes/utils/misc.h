/*
 * Copyright by Paul Rudolph
 * Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
 * https://www.leibniz-hki.de/en/applied-systems-biology.html
 * HKI-Center for Systems Biology of Infection
 * Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
 * Adolf-Reichwein-Straße 23, 07745 Jena, Germany
 *
 * This code is licensed under BSD 2-Clause
 * See the LICENSE file provided with this code for the full license.
 */

#ifndef SBM_MISC_H
#define SBM_MISC_H

#include <iomanip>
#include <map>
#include <sstream>
#include <vector>

namespace sbm::util {


template<typename T>
class VvIterator : public std::iterator<std::bidirectional_iterator_tag, T>{
  public:

    static VvIterator<T> begin(std::vector<std::vector<T>>& vv) {
        return VvIterator(&vv, 0, 0);
    }
    static VvIterator<T> end(std::vector<std::vector<T>>& vv) {
        return VvIterator(&vv, vv.size(), 0);
    }

    VvIterator() = default;
    // ++prefix operator
    VvIterator& operator++()
    {
        // If we haven't reached the end of this sub-vector.
        if (idx_inner_ + 1 < (*vv_)[idx_outer_].size())
        {
            // Go to the next element.
            ++idx_inner_;
        }
        else
        {
            // Otherwise skip to the next sub-vector, and keep skipping over empty
            // ones until we reach a non-empty one or the end.
            do
            {
                ++idx_outer_;
            } while (idx_outer_ < (*vv_).size() && (*vv_)[idx_outer_].empty());

            // Go to the start of this vector.
            idx_inner_ = 0;
        }
        return *this;
    }
    // --prefix operator
    VvIterator& operator--()
    {
        // If we haven't reached the start of this sub-vector.
        if (idx_inner_ > 0)
        {
            // Go to the previous element.
            --idx_inner_;
        }
        else
        {
            // Otherwise skip to the previous sub-vector, and keep skipping over empty
            // ones until we reach a non-empty one.
            do
            {
                --idx_outer_;
            } while ((*vv_)[idx_outer_].empty());

            // Go to the end of this vector.
            idx_inner_ = (*vv_)[idx_outer_].size() - 1;
        }
        return *this;
    }
    // postfix++ operator
    VvIterator operator++(int)
    {
        T retval = *this;
        ++(*this);
        return retval;
    }
    // postfix-- operator
    VvIterator operator--(int)
    {
        T retval = *this;
        --(*this);
        return retval;
    }
    bool operator==(const VvIterator& other) const
    {
        return other.vv_ == vv_ && other.idx_outer_ == idx_outer_ && other.idx_inner_ == idx_inner_;
    }
    bool operator!=(const VvIterator &other) const
    {
        return *this != other;
    }
    const T& operator*() const
    {
        return *this;
    }
    T& operator*()
    {
        return (*vv_)[idx_outer_][idx_inner_];
    }
    const T& operator->() const
    {
        return *this;
    }
    T& operator->()
    {
        return *this;
    }

  private:
    VvIterator(std::vector<std::vector<T>>* vv,
                std::size_t idx_outer,
                std::size_t idx_inner)
        : vv_(vv), idx_outer_(idx_outer), idx_inner_(idx_inner) {}

    std::vector<std::vector<T>>* vv_ = nullptr;
    std::size_t idx_outer_ = 0;
    std::size_t idx_inner_ = 0;
};



} // namespace sbm::util

#endif //SBM_MISC_H
