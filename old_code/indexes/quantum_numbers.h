#ifndef JULY_QUANTUM_NUMBERS_H
#define JULY_QUANTUM_NUMBERS_H

// Quantum_Numbers is vector with size N + (N - 1).
// First N values is particles' spins, other (N - 1) values is the sums of spins,
// produced during additions of particles' spins.
// Order of summation is written in addition_scheme vector of Indexes class.

struct Quantum_Numbers {
public:
    unsigned int representation;

    const int & operator()(unsigned long i) const {
        if (i < mults.size()) {
            return mults[i];
        } else {
            return sum_mults[i - mults.size()];
        }
    }

    // NB: decreasing order of back_mult!
    bool operator<(const Quantum_Numbers &other) const {
        return (this->representation < other.representation) ||
               ((this->representation == other.representation) && (this->back_mult() > other.back_mult()));
    }

    Quantum_Numbers(const std::vector<int> & mults_, unsigned int repr) : mults(mults_), representation(repr) {
        sum_mults.reserve(mults.size() - 1);
    }

    Quantum_Numbers& operator=(Quantum_Numbers&& other)  noexcept {
        representation = other.representation;
        sum_mults = std::move(other.sum_mults);
        return *this;
    }

    Quantum_Numbers(const Quantum_Numbers&) = default;
    Quantum_Numbers& operator= (const Quantum_Numbers&) = delete;
    Quantum_Numbers(Quantum_Numbers&&) noexcept = default;
    ~Quantum_Numbers() noexcept = default;

    bool empty() const {
        return sum_mults.empty();
    }

    unsigned long size() const {
        return sum_mults.size();
    }

    int back_mult() const {
        if (this->empty()) {
            return 1;
        } else {
            return sum_mults.back();
        }
    }

    void push_back(int m) {
        sum_mults.push_back(m);
    }

    void print() {
        for (auto m : sum_mults) {
            std::cout << m << " ";
        }
        std::cout << std::endl;
    }

private:
    std::vector<int> sum_mults;
    const std::vector<int> & mults;
};


#endif //JULY_QUANTUM_NUMBERS_H
