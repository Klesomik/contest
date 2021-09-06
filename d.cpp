#include <bits/stdc++.h>

// https://github.com/dmkz/competitive-programming/blob/master/e-olymp.com/0317.cpp
const long double PI = std::acos(-1.0L);

struct UInt {
    static const int BASE = (int)1e9; // Основание системы счисления
    static const int WIDTH = 9;       // Количество десятичных цифр, которые хранятся в одной цифре
    
    // Вектор под цифры числа:
    std::vector<int> digits; 
    
    // Конструкторы
    UInt(int64_t number = 0);
    UInt(const std::string& s);
    UInt(const std::vector<int>& digits);
    
    // Методы нормализации и сравнения:
    UInt& normalize(); // удаление лидирующих нулей и проверка на принадлежность цифр диапазону [0, BASE)
    int compare(const UInt& other) const; // Сравнение (меньше = -1, равно = 0, больше = 1)
    
    // Методы умножения:
    UInt slow_mult(const UInt& other) const; // Медленное произведение (работает довольно быстро на числах небольшой длины)
    UInt fast_mult(const UInt& other) const; // Быстрое произведение (на основе Быстрого Преобразования Фурье комплексные числа)
    UInt mult(const UInt& other) const; // Комбинированный метод умножения на основе экспериментальных данных
    
    // Метод деления:
    std::pair<UInt, UInt> div_mod(const UInt& other) const; // Целая часть и остаток от деления
    
    // Операторы:
    UInt& operator+=(const uint32_t num);     // Прибавление короткого
    UInt& operator+=(const UInt& other); // Прибавление длинного
    UInt& operator-=(const int num);     // Вычитание короткого
    UInt& operator-=(const UInt& other); // Вычитание длинного
    UInt& operator*=(const uint32_t num);     // Умножение на короткое
    UInt& operator*=(const UInt& other); // Умножение на длинное
    UInt& operator/=(const uint32_t num);     // Деление на короткое
    UInt& operator/=(const UInt& other); // Деление на длинное
    UInt& operator%=(const UInt& other); // Остаток от деления на длинное

    uint32_t to_uint32t() const
    {
        uint32_t result = 0;
        for (auto it = digits.rbegin(); it != digits.rend(); it++)
        {
            result *= BASE;
            result += (*it);
        }
        return result;
    }
};

std::istream& operator>>(std::istream&, UInt&); // Ввод из потока
std::ostream& operator<<(std::ostream&, const UInt&); // Вывод в поток

UInt pow(UInt, int); // Возведение в степень
UInt gcd(UInt, UInt); // Наибольший общий делитель

UInt operator+(const UInt&, const UInt&);
UInt operator-(const UInt&, const UInt&);
UInt operator*(const UInt&, const UInt&);
UInt operator/(const UInt&, const UInt&);
UInt operator%(const UInt&, const UInt&);

UInt operator+(const UInt&, const int);
UInt operator+(const int, const UInt&);
UInt operator-(const UInt&, const int);
UInt operator*(const UInt&, const int);
UInt operator*(const int, const UInt&);
UInt operator/(const UInt&, const int);
UInt operator^(const UInt&, const int); // возведение в степень

bool operator<(const UInt&, const UInt&);
bool operator>(const UInt&, const UInt&);
bool operator<=(const UInt&, const UInt&);
bool operator>=(const UInt&, const UInt&);
bool operator==(const UInt&, const UInt&);
bool operator!=(const UInt&, const UInt&);

UInt& UInt::normalize() {
    while (digits.back() == 0 && (int)digits.size() > 1) digits.pop_back();
    for (auto d : digits) assert(0 <= d && d < BASE);
    return *this;
}   

// Конструктор от короткого целого
UInt::UInt(int64_t number) {
    assert(number >= 0);
    do {
        digits.push_back(number % BASE);
        number /= BASE;
    } while (number > 0);
    normalize();
}

// Конструктор от вектора из цифр:
UInt::UInt(const std::vector<int>& digits) : digits(digits) { 
    normalize();
}

// Конструктор от строчки:
UInt::UInt(const std::string& s) {
    const int size = (int)s.size();
    for (int idGroup = 1, nGroups = size / WIDTH; idGroup <= nGroups; ++idGroup) {            
        digits.push_back(std::stoi(s.substr(size-idGroup * WIDTH, WIDTH)));
    }
    if (size % WIDTH != 0) {
        digits.push_back(std::stoi(s.substr(0, size % WIDTH)));
    }
    if (digits.empty()) {
        digits.push_back(0);
    }
    normalize();
}

// Прибавление короткого:
UInt& UInt::operator+=(const uint32_t num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this += UInt(num);
    } 
    uint32_t rem = num;
    for (int i = 0; rem > 0; ++i) {
        if (i >= (int)digits.size()) digits.push_back(0);
        rem += digits[i];
        if (rem >= BASE) {
            digits[i] = rem - BASE;
            rem = 1;
        } else {
            digits[i] = rem;
            rem = 0;
        }
    }
    return this->normalize();
}

// Прибавление длинного:
UInt& UInt::operator+=(const UInt& other) {
    if (other.digits.size() == 1u) {
        return *this += other.digits[0];
    }
    const int s1 = this->digits.size();
    const int s2 = other.digits.size();
    int rem = 0;
    for (int i = 0; i < s1 || i < s2 || rem > 0; ++i) {
        int d1 = i < s1 ? this->digits[i] : (digits.push_back(0), 0);
        int d2 = i < s2 ? other.digits[i] : 0;
        rem += d1 + d2;
        auto div = rem / BASE;
        digits[i] = rem - div * BASE;
        rem = div;
    }
    return this->normalize();
}

// Вычитание короткого:
UInt& UInt::operator-=(const int num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this -= UInt(num);
    }
    int rem = -num;
    for (int i = 0; i < (int)digits.size() && rem < 0; ++i) {
        rem += digits[i];
        if (rem < 0) { // Занимаем разряд
            digits[i] = rem + BASE;
            rem = -1;
        } else {
            digits[i] = rem;
            rem = 0;
        }
    }
    assert(rem == 0);
    return this->normalize();
}

// Вычитание длинного:
UInt& UInt::operator-=(const UInt& other) {
    if (other.digits.size() == 1u) {
        return *this -= other.digits[0];
    }
    const int s1 = this->digits.size();
    const int s2 = other.digits.size();
    assert(s1 >= s2);
    int rem = 0;
    for (int i = 0; i < s1; ++i) {
        int d2 = i < s2 ? other.digits[i] : 0;
        rem += this->digits[i] - d2;
        if (rem < 0) {
            digits[i] = rem + BASE;
            rem = -1;
        } else {
            digits[i] = rem;
            rem = 0;
            if (i >= s2) break;
        }
    }
    assert(rem == 0); // Иначе *this < other
    return this->normalize();
}

// Умножение на короткое:
UInt& UInt::operator*=(const uint32_t num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this *= UInt(num);
    }
    int64_t rem = 0;
    for (auto& d : digits) {
        rem += 1LL * d * num;
        auto div = rem / BASE;
        d = rem - div * BASE;
        rem = div;
    }
    if (rem > 0) digits.push_back(rem);
    return this->normalize();
}

// Медленное произведение:
UInt UInt::slow_mult(const UInt& other) const {
    if (other.digits.size() == 1u) {
        return *this * other.digits[0];
    }
    const int s1 = (int)this->digits.size();
    const int s2 = (int)other.digits.size();
    std::vector<int> temp(s1+s2);
    for (int i = 0; i < s1; ++i) {
        int64_t rem = 0;
        for (int j = 0; j < s2; ++j) {
            rem += temp[i+j] + 1LL * this->digits[i] * other.digits[j];
            auto div = rem / BASE;
            temp[i+j] = rem - div * BASE;
            rem = div;
        }
        if (rem > 0) temp[i+s2] += rem;
    }
    return UInt(temp);
}

// Быстрое умножение на основе быстрого преобразования Фурье:
UInt UInt::fast_mult(const UInt& other) const {
    if (other.digits.size() == 1u) {
        return *this * other.digits[0];
    }
    
    // Разворот битов в числе num:
    std::function<int(int, int)> reverse = [](int number, int nBits) {
        int res = 0;
        for (int i = 0; i < nBits; ++i) {
            if (number & (1 << i)) {
                res |= 1 << (nBits-1-i);
            }
        }
        return res;
    };
    
    typedef std::complex<long double> complex;
    // Быстрое преобразование Фурье:
    std::function<void(std::vector<complex>&, bool)> fft = [&reverse](std::vector<complex> & a, bool invert) {
        const int n = (int)a.size();
        int nBits = 0;
        while ((1 << nBits) < n) ++nBits;

        for (int i = 0; i < n; ++i) {
            if (i < reverse(i, nBits)) {
                std::swap(a[i], a[reverse(i, nBits)]);
            }
        }

        for (int len = 2; len <= n; len <<= 1) {
            auto ang = 2*PI / len * (invert ? -1 : 1);
            complex wlen (std::cos(ang), std::sin(ang));
            for (int i = 0; i < n; i += len) {
                complex w(1);
                for (int j = 0; j < len / 2; ++j) {
                    complex u = a[i+j];
                    complex v = a[i+j+len / 2] * w;
                    a[i+j] = u + v;
                    a[i+j+len/2] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert) {
            for (int i = 0; i < n; ++i) {
                a[i] /= n;
            }
        }
    };
    
    // Подготавливаем вектора из комплексных коэффициентов fa и fb:
    // Так как происходит потеря точности из-за арифметики с плавающей точкой, основание системы необходимо понизить:
    assert(BASE == 1000 * 1000 * 1000);
    std::function<std::vector<complex>(const UInt&)> prepare = [](const UInt& number) {
        std::vector<complex> result; 
        result.reserve(3 * number.digits.size());
        for (auto d : number.digits) {
            result.push_back(d % 1000);
            result.push_back(d / 1000 % 1000);
            result.push_back(d / 1000000);    
        }
        return result;
    };
    
    auto fa = prepare(*this);
    auto fb = prepare(other);
    
    // Округляем размер векторов до ближайшей степени двойки:
    int n = 1;
    while (n < (int)std::max(fa.size(), fb.size())) n *= 2;
    n *= 2;
    fa.resize(n);
    fb.resize(n);
    
    // Вызываем прямое преобразование Фурье:
    fft (fa, false);
    fft (fb, false);
    // Перемножаем результаты:
    for (int i = 0; i < n; ++i) {
        fa[i] *= fb[i];
    }
    // Вызываем обратное преобразование Фурье:
    fft (fa, true);
    // Копируем ответ с округлениями:
    std::vector<int64_t> temp(n);
    for (int i = 0; i < (int)fa.size(); ++i) {
        temp[i] = int64_t (fa[i].real() + 0.5);
    }
    // Не забываем про переносы в старшие разряды:
    int64_t carry = 0;
    for (int i = 0; i < n || carry > 0; ++i) {
        if (i >= n) temp.push_back(0);
        temp[i] += carry;
        carry = temp[i] / 1000;
        temp[i] -= carry * 1000;
        assert(temp[i] >= 0);
    }
    // Формируем ответ:
    std::vector<int> res;
    res.reserve(this->digits.size() + other.digits.size());
    
    for (int i = 0; i < n; i += 3) {
        int c = temp[i];
        int b = i+1 < n ? temp[i+1] : 0;
        int a = i+2 < n ? temp[i+2] : 0;
        res.push_back(c + 1000 * (b + 1000 * a));
    }
    return UInt(res);
}

// Комбинированный метод умножения:
UInt UInt::mult(const UInt& other) const {
// Выбор метода умножения:
    int len1 = (int)this->digits.size();
    int len2 = (int)other.digits.size();
    int temp = 3 * std::max(len1, len2);
    int pow = 1;
    while (pow < temp) pow *= 2;
    pow *= 2;
    int op1 = len1 * len2;
    int op2 = 3 * pow * std::log(pow) / std::log(2);
    return op1 >= 15 * op2 ? fast_mult(other) : slow_mult(other);
}

// Деление на короткое:
UInt& UInt::operator/=(const uint32_t num) {
    assert(num > 0);
    if (num >= BASE) {
        return *this /= UInt(num);
    }
    int64_t rem = 0;
    for (int j = (int)digits.size()-1; j >= 0; --j) {
        (rem *= BASE) += digits[j];
        auto div = rem / num;
        digits[j] = div;
        rem -= div * num;
    }
    return this->normalize();
}

// Остаток от деления на короткое:
uint32_t operator%(const UInt& a, const uint32_t num) {
    assert(num > 0);
    int64_t rem = 0;
    for (int i = (int)a.digits.size()-1; i >= 0; --i) {
        ((rem *= UInt::BASE) += a.digits[i]) %= num;
    }
    return rem;
}

// Целая часть и остаток от деления:
std::pair<UInt, UInt> UInt::div_mod(const UInt& other) const {
    if (other.digits.size() == 1u) {
        return {std::move(*this / other.digits[0]), *this % other.digits[0]};
    }
    const int norm = BASE / (other.digits.back() + 1);
    const UInt a = *this * norm;
    const UInt b = other * norm;
    const int a_size = (int)a.digits.size();
    const int b_size = (int)b.digits.size();
    UInt q, r;
    q.digits.resize(a_size);
    for (int i = a_size - 1; i >= 0; --i) {
        r *= BASE;
        r += a.digits[i];
        int s1 = (int)r.digits.size() <= b_size ? 0 : r.digits[b_size];
        int s2 = (int)r.digits.size() <= b_size - 1 ? 0 : r.digits[b_size - 1];
        int d = (1LL * BASE * s1 + s2) / b.digits.back();
        auto temp = b * d;
        while (r < temp) {
            r += b;
            --d;
        }
        r -= temp;
        q.digits[i] = d;
    }
    return {std::move(q.normalize()), std::move(r /= norm)};
}

// Сравнение: result < 0 (меньше), result == 0 (равно), result > 0 (больше)
int UInt::compare(const UInt& other) const {
    if (this->digits.size() > other.digits.size()) return 1;
    if (this->digits.size() < other.digits.size()) return -1;
    for (int i = (int)digits.size()-1; i >= 0; --i) {
        if (this->digits[i] > other.digits[i]) return 1;
        if (this->digits[i] < other.digits[i]) return -1;
    }
    return 0;
}

// Операторы сравнения:
bool operator< (const UInt& a, const UInt& b) { return a.compare(b) < 0; }
bool operator> (const UInt& a, const UInt& b) { return a.compare(b) > 0; }
bool operator==(const UInt& a, const UInt& b) { return a.compare(b) == 0; }
bool operator<=(const UInt& a, const UInt& b) { return a.compare(b) <= 0; }
bool operator>=(const UInt& a, const UInt& b) { return a.compare(b) >= 0; }
bool operator!=(const UInt& a, const UInt& b) { return a.compare(b) != 0; }

// Ввод из потока:
std::istream& operator>>(std::istream& is, UInt& number) {
    std::string s;
    is >> s;
    number = UInt(s);
    return is;
}

// Вывод в поток:
std::ostream& operator<<(std::ostream& os, const UInt& number) {
    os << number.digits.back();
    for (int i = (int)number.digits.size()-2; i >= 0; --i) {
        os << std::setw(UInt::WIDTH) << std::setfill('0') << number.digits[i];
    }
    return os << std::setfill(' ');
}

// Сумма:
UInt operator+(const UInt& a, const UInt& b) { 
    return UInt(a) += b; 
}

// Разность:
UInt operator-(const UInt& a, const UInt& b) { 
    return UInt(a) -= b; 
}

// Произведение:
UInt operator*(const UInt& a, const UInt& b) { 
    return a.mult(b);
}

// Деление:
UInt operator/(const UInt& a, const UInt& b) {
    return a.div_mod(b).first;
}

// Взятие остатка:
UInt operator%(const UInt& a, const UInt& b) {
    return a.div_mod(b).second;
}

// Умножение:
UInt& UInt::operator*=(const UInt& other) {
    return *this = *this * other;
}

// Деление с присваиванием:
UInt& UInt::operator/=(const UInt& other) {
    return *this = *this / other;
}

// Взятие остатка с присваиванием:
UInt& UInt::operator%=(const UInt& other) {
    return *this = *this % other;
}

UInt operator+(const UInt& a, const int b) { return UInt(a) += b; }
UInt operator+(const int a, const UInt& b) { return b * a; }
UInt operator-(const UInt& a, const int b) { return UInt(a) -= b; }
UInt operator*(const UInt& a, const int b) { return UInt(a) *= b; }
UInt operator*(const int a, const UInt& b) { return b * a; }
UInt operator/(const UInt& a, const int b) { return UInt(a) /= b; }
UInt operator^(const UInt& a, const int n) { return pow(a, n); } // Возведение в степень

// Возведение в степень:
UInt pow(UInt a, int n) {
    UInt res = 1;
    while (n > 0) {
        if (n % 2 != 0) res *= a;
        a *= a;
        n /= 2;
    }
    return res;
}

// Наибольший общий делитель:
UInt gcd(UInt a, UInt b) {
    while (b != 0) {
        auto rem = a % b;
        a = b;
        b = rem;
    }
    return a;
}

uint64_t mod(uint64_t a, uint64_t b)
{
    return (a % b + b) % b;
}

int64_t mod_inv(int64_t a, int64_t b)
{
    return (a % b + b) % b;
}

uint32_t p = 0;

/**
 * Field Zp, p -- prime
 */
struct Zp
{
    uint32_t value;

    Zp(int64_t value_ = 0): value(mod_inv(value_, p))
    {
    }

    Zp& operator*=(const Zp& other)
    {
        value = mod(1ULL * value * other.value, p);
        return (*this);
    }

    Zp& operator+=(const Zp& other)
    {
        value = mod(1ULL * value + other.value, p);
        return (*this);
    }

    Zp& operator-=(const Zp& other)
    {
        int64_t tmp = other.value;
        tmp = -tmp;
        return (*this) += mod_inv(tmp, p);
    }

    uint32_t get_ord() const
    {
        return p - 1;
    }

    Zp get_e() const
    {
        return Zp(1);
    }
};

std::ostream& operator<<(std::ostream& os, const Zp& zp)
{
    return os << zp.value;
}

std::istream& operator>>(std::istream& is, Zp& zp)
{
    return is >> zp.value;
}

bool operator==(const Zp& lhs, const Zp& rhs)
{
    return lhs.value == rhs.value;
}

Zp operator*(const Zp& lhs, const Zp& rhs)
{
    Zp result = lhs;
    result *= rhs;
    return result;
}

Zp pow(const Zp& a, uint32_t n)
{
    if (n == 0)
    {
        return a.get_e();
    }
    if (n % 2 == 0)
    {
        Zp tmp = pow(a, n / 2);
        return tmp * tmp;
    }
    return a * pow(a, n - 1);
}

Zp get_inv(const Zp& a)
{
    return pow(a, p - 2);
}

Zp operator/(const Zp& lhs, const Zp& rhs)
{
    return lhs * pow(rhs, p - 2);
}

std::istream& operator>>(std::istream& is, std::vector<Zp>& buffer)
{
    std::string line;
    std::getline(is, line, '\n');
    if (line == "")
    {
        std::getline(is, line, '\n');
    }
    std::istringstream ss(line);
    for (int64_t coeff = 0; ss >> coeff;)
    {
        buffer.push_back(Zp(coeff));
    }
    return is;
}

std::ostream& operator<<(std::ostream& os, const std::vector<Zp>& buffer)
{
    for (const auto &it : buffer)
    {
        os << it << ' ';
    }
    return os;
}

std::vector<Zp> f;
uint32_t n = 0;

/**
 * Field Fq = Zp[x]/(f), q = p^n
 */
struct Fq
{
    std::vector<Zp> coeffs;

    Fq(): coeffs()
    {
    }

    Fq(const std::vector<Zp>& coeffs_): coeffs(coeffs_)
    {
    }

    Fq& operator*=(const Fq& other)
    {
        assert(coeffs.size() == n);
        assert(other.coeffs.size() == n);
        std::vector<Zp> buffer(n * 2 - 1, Zp(0));
        for (uint32_t i = 0; i < n; i++)
        {
            for (uint32_t j = 0; j < n; j++)
            {
                buffer[i + j] += coeffs[i] * other.coeffs[j];
            }
        }
        while (buffer.size() >= f.size())
        {
            Zp alpha = buffer.back() / f.back();
            for (uint32_t i = 0; i < f.size(); i++)
            {
                buffer[buffer.size() - 1 - i] -= alpha * f[f.size() - 1 - i];
            }
            buffer.pop_back();
        }
        coeffs = buffer;
        return (*this);
    }

    uint32_t get_ord() const
    {
        return p * n;
    }

    Fq get_e() const
    {
        Fq tmp(std::vector<Zp>(n, Zp(0)));
        tmp.coeffs[0] = Zp(1);
        return tmp;
    }
};

std::ostream& operator<<(std::ostream& os, const Fq& fq)
{
    return os << fq.coeffs;
}

std::istream& operator>>(std::istream& is, Fq& fq)
{
    is >> fq.coeffs;
    fq.coeffs.resize(n);
    return is;
}

bool operator==(const Fq& lhs, const Fq& rhs)
{
    return lhs.coeffs == rhs.coeffs;
}

Fq operator*(const Fq& lhs, const Fq& rhs)
{
    Fq result = lhs;
    result *= rhs;
    return result;
}

Fq pow(const Fq& a, uint32_t n)
{
    if (n == 0)
    {
        return a.get_e();
    }
    if (n % 2 == 0)
    {
        Fq tmp = pow(a, n / 2);
        return tmp * tmp;
    }
    return a * pow(a, n - 1);
}

namespace FqTools
{
    void normalize(Fq& a)
    {
        while ((a.coeffs.size() > 1) && (a.coeffs.back() == Zp(0)))
        {
            a.coeffs.pop_back();
        }
    }

    bool is_zero(const Fq& a)
    {
        return (a.coeffs.size() == 1) && (a.coeffs[0] == Zp(0));
    }

    Fq add(Fq a, Fq b)
    {
        if (a.coeffs.size() < b.coeffs.size())
        {
            a.coeffs.resize(b.coeffs.size());
        }
        for (uint32_t i = 0; i < b.coeffs.size(); i++)
        {
            a.coeffs[i] += b.coeffs[i];
        }
        normalize(a);
        return a;
    }

    Fq sub(Fq a, Fq b)
    {
        if (a.coeffs.size() < b.coeffs.size())
        {
            a.coeffs.resize(b.coeffs.size());
        }
        for (uint32_t i = 0; i < b.coeffs.size(); i++)
        {
            a.coeffs[i] -= b.coeffs[i];
        }
        normalize(a);
        return a;
    }

    Fq mul(Fq a, Fq b)
    {
        std::vector<Zp> buffer(a.coeffs.size() + b.coeffs.size() - 1, Zp(0));
        for (uint32_t i = 0; i < a.coeffs.size(); i++)
        {
            for (uint32_t j = 0; j < b.coeffs.size(); j++)
            {
                buffer[i + j] += a.coeffs[i] * b.coeffs[j];
            }
        }
        Fq c(buffer);
        normalize(c);
        return c;
    }

    Fq div(Fq a, Fq b)
    {
        normalize(a);
        normalize(b);
        std::vector<Zp> result(a.coeffs.size() - b.coeffs.size() + 1);
        while (a.coeffs.size() >= b.coeffs.size())
        {
            Zp alpha = a.coeffs.back() / b.coeffs.back();
            result[a.coeffs.size() - b.coeffs.size()] = alpha;
            for (uint32_t i = 0; i < b.coeffs.size(); i++)
            {
                a.coeffs[a.coeffs.size() - 1 - i] -= alpha * b.coeffs[b.coeffs.size() - 1 - i];
            }
            a.coeffs.pop_back();
        }
        Fq c(result);
        normalize(c);
        return c;
    }

    uint32_t deg(Fq a)
    {
        for (uint32_t i = 0; i < a.coeffs.size(); i++)
        {
            if (!(a.coeffs[a.coeffs.size() - 1 - i] == Zp(0)))
            {
                return a.coeffs.size() - 1 - i;
            }
        }
        assert(false);
    }

    Fq gcdex(Fq a, Fq b, Fq& u, Fq& v)
    {
        normalize(a);
        normalize(b);
        if (is_zero(a))
        {
            v.coeffs.push_back(Zp(1));
            return b;
        }
        if (is_zero(b))
        {
            u.coeffs.push_back(Zp(1));
            return a;
        }
        if (deg(b) >= deg(a))
        {
            Fq q = div(b, a);
            Fq r = sub(b, mul(a, q));
            assert(b == add(mul(a, q), r));
            Fq h = gcdex(a, r, u, v);
            u = sub(u, mul(q, v));
            return h;
        }
        return gcdex(b, a, v, u);
    }

    Fq get_inv(Fq a)
    {
        Fq u, v;
        Fq h = gcdex(a, Fq(f), u, v);
        normalize(h);
        normalize(u);
        normalize(v);
        assert(is_zero(sub(a, mul(h, div(a, h)))));
        assert(add(mul(a, u), mul(Fq(f), v)) == h);
        if (!(h.coeffs[0] == Zp(1)))
        {
            for (auto &it : u.coeffs)
            {
                it = it / h.coeffs[0];
            }
            for (auto &it : v.coeffs)
            {
                it = it / h.coeffs[0];
            }
            h.coeffs[0] = Zp(1);
        }
        assert((h.coeffs.size() == 1) && (h.coeffs[0] == Zp(1)));
        assert(is_zero(sub(a, mul(h, div(a, h)))));
        assert(add(mul(a, u), mul(Fq(f), v)) == h);
        u.coeffs.resize(n);
        return u;
    }
}

char number_to_char(int number)
{
    if (number == 62)
        return ' ';
    if (number == 63)
        return '.';
    if (number >= 0 && number <= 9)
        return '0' + number;
    number -= '9' - '0' + 1;
    if (number >= 0 && number <= 25)
        return 'A' + number;
    number -= 'Z' - 'A' + 1;
    if (number >= 0 && number <= 25)
        return 'a' + number;
    assert(false);
}

std::vector<Fq> decrypt_el_gamal(uint32_t a, const std::vector<std::pair<Fq, Fq>>& buffer)
{
    std::vector<Fq> result;
    for (const auto &it : buffer)
    {
        Fq r_i = it.first;
        Fq m_i = it.second;

        Fq k = pow(r_i, a);
        Fq k_inv = FqTools::get_inv(k);
        assert(k_inv.coeffs.size() == n);
        assert(k * k_inv == k.get_e());

        result.push_back(m_i * k_inv);
    }

    return result;
}

int main()
{
    /*std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.tie(0);*/

    /* Init global vars */
    std::cin >> p;
    std::cin >> f;
    n = f.size() - 1;

    /* a -- private key */
    uint32_t a = 0;
    std::cin >> a;

    std::vector<std::pair<Fq, Fq>> buffer;
    while (true)
    {
        Fq r_i;
        if (!(std::cin >> r_i))
        {
            break;
        }

        Fq m_i;
        std::cin >> m_i;

        buffer.push_back({ r_i, m_i });
    }

    std::vector<Fq> result = decrypt_el_gamal(a, buffer);

    UInt number = 0;
    UInt base = 1;
    for (const auto &fq : result)
    {
        for (const auto &it : fq.coeffs)
        {
            number += base * UInt(it.value);
            base *= UInt(p);
        }
    }

    std::string message;
    while (number > 0)
    {
        message.push_back(number_to_char(number % 64));
        number /= 64;
    }

    std::cout << message << '\n';

    return 0;
}
