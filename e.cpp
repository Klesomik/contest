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
    UInt& operator+=(const int num);     // Прибавление короткого
    UInt& operator+=(const UInt& other); // Прибавление длинного
    UInt& operator-=(const int num);     // Вычитание короткого
    UInt& operator-=(const UInt& other); // Вычитание длинного
    UInt& operator*=(const int num);     // Умножение на короткое
    UInt& operator*=(const UInt& other); // Умножение на длинное
    UInt& operator/=(const uint32_t num);     // Деление на короткое
    UInt& operator/=(const UInt& other); // Деление на длинное
    UInt& operator%=(const UInt& other); // Остаток от деления на длинное
};

std::istream& operator>>(std::istream&, UInt&); // Ввод из потока
std::ostream& operator<<(std::ostream&, const UInt&); // Вывод в поток

UInt pow(UInt, UInt); // Возведение в степень

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
UInt operator^(const UInt&, const UInt); // возведение в степень

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
UInt& UInt::operator+=(const int num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this += UInt(num);
    } 
    int rem = num;
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
UInt& UInt::operator*=(const int num) {
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
UInt operator^(const UInt& a, const UInt n) { return pow(a, n); } // Возведение в степень

// Возведение в степень:
UInt pow(UInt a, UInt n) {
    UInt res = 1;
    while (n > 0) {
        if (n % 2 != 0) res *= a;
        a *= a;
        n /= 2;
    }
    return res;
}

UInt p = (UInt(2) ^ UInt(256)) - (UInt(2) ^ UInt(224)) + (UInt(2) ^ UInt(192)) + (UInt(2) ^ UInt(96)) - 1;

/**
 * Field Zp, p -- prime
 */
struct Zp
{
    UInt value;

    Zp(UInt value_ = 0): value(value_ % p)
    {
    }

    Zp& operator*=(const Zp& other)
    {
        value = (value * other.value) % p;
        return (*this);
    }

    Zp& operator+=(const Zp& other)
    {
        value = (value + other.value) % p;
        return (*this);
    }

    Zp& operator-=(const Zp& other)
    {
        if (value < other.value)
        {
            assert((p + value) >= other.value);
            value = p + value - other.value;
            assert((value >= UInt(0)) && (value < p));
        }
        else
        {
            value = value - other.value;
            assert((value >= UInt(0)) && (value < p));
        }
        return (*this);
    }

    Zp pow(UInt n) const
    {
        if (n == 0)
        {
            return Zp(1);
        }
        if (n % 2 == 0)
        {
            Zp tmp = pow(n / 2);
            tmp *= tmp;
            return tmp;
        }
        Zp tmp = (*this);
        tmp *= pow(n - 1);
        return tmp;
    }

    Zp get_inv() const
    {
        return pow(p - 2);
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

bool operator!=(const Zp& lhs, const Zp& rhs)
{
    return !(lhs == rhs);
}

Zp operator+(const Zp& lhs, const Zp& rhs)
{
    Zp result = lhs;
    result += rhs;
    return result;
}

Zp operator-(const Zp& lhs, const Zp& rhs)
{
    Zp result = lhs;
    result -= rhs;
    return result;
}

Zp operator*(const Zp& lhs, const Zp& rhs)
{
    Zp result = lhs;
    result *= rhs;
    return result;
}

Zp operator/(const Zp& lhs, const Zp& rhs)
{
    return lhs * rhs.get_inv();
}

/* NIST P-256 */
UInt a("3"); // -3
UInt b("41058363725152142129326129780047268409114441015993725554835256314039467401291");

struct EZp
{
    Zp x, y;
    bool is_inf;

    EZp(): x(), y(), is_inf(false)
    {
    }

    EZp(const Zp& x_): x(x_), y(), is_inf(false)
    {
        set_y();
    }

    EZp(const Zp& x_, const Zp& y_): x(x_), y(y_), is_inf(false)
    {
        assert((y * y) == (x * x * x - (Zp(a) * x) + Zp(b)));
    }

    void set_y()
    {
        UInt k = (p + UInt(1)) / UInt(4);
        Zp result = x.pow(3) - (Zp(a) * x) + Zp(b);
        y = result.pow(k);
        assert(!is_inf);
        assert((y * y) == (x * x * x - (Zp(a) * x) + Zp(b)));
    }
};

EZp get_e()
{
    EZp tmp;
    tmp.is_inf = true;
    return tmp;
}

std::ostream& operator<<(std::ostream& os, const EZp& ezp)
{
    if (ezp.is_inf)
    {
        return os << 'Z';
    }
    assert((ezp.y * ezp.y) == (ezp.x * ezp.x * ezp.x - (Zp(a) * ezp.x) + Zp(b)));
    return os << ezp.x << ' ' << ezp.y;
}

std::istream& operator>>(std::istream& is, EZp& ezp)
{
    is >> ezp.x >> ezp.y;
    return is;
}

EZp operator+(const EZp& lhs, const EZp& rhs)
{
    if (lhs.is_inf)
    {
        return rhs;
    }
    if (rhs.is_inf)
    {
        return lhs;
    }
    if (lhs.x != rhs.x)
    {
        Zp k = (rhs.y - lhs.y) / (rhs.x - lhs.x);
        Zp x = k * k - lhs.x - rhs.x;
        Zp y = k * (lhs.x - x) - lhs.y;
        return EZp(x, y);
    }
    else if (lhs.y == rhs.y)
    {
        if (lhs.y == UInt(0))
        {
            return get_e();
        }
        Zp k = (lhs.x * lhs.x * Zp(3) - Zp(a)) / (lhs.y * Zp(2));
        Zp x = k * k - Zp(2) * lhs.x;
        Zp y = k * (lhs.x - x) - lhs.y;
        return EZp(x, y);
    }
    else
    {
        return get_e();
    }
}

EZp operator*(const EZp& lhs, UInt n)
{
    if (n == 0)
    {
        return get_e();
    }
    if (n % 2 == 0)
    {
        EZp tmp = lhs * (n / 2);
        return tmp + tmp;
    }
    return lhs + (lhs * (n - 1));
}

bool operator==(const EZp& lhs, const EZp& rhs)
{
    if (lhs.is_inf != rhs.is_inf)
    {
        return false;
    }
    if (lhs.is_inf == rhs.is_inf)
    {
        return true;
    }
    return (lhs.x == rhs.x) && (lhs.y == rhs.y);
}

EZp g(UInt("48439561293906451759052585252797914202762949526041747995844080717082404635286"), UInt("36134250956749795798585127919587881956611106672985015071877198253568414405109"));
UInt ord("115792089210356248762697446949407573529996955224135760342422259061068512044369");

int char_to_number(char symbol)
{
    if (symbol >= '0' && symbol <= '9')
        return symbol - '0';
    if (symbol >= 'A' && symbol <= 'Z')
        return symbol - 'A' + ('9' - '0' + 1);
    if (symbol >= 'a' && symbol <= 'z')
        return symbol - 'a' + ('9' - '0' + 1) + ('Z' - 'A' + 1);
    if (symbol == '_')
        return 62;
    if (symbol == '.')
        return 63;
    assert(false);
}

std::vector<EZp> scan_message()
{
    int n = 0;
    std::cin >> n;

    std::vector<EZp> message;
    for (int i = 0; i < n; i++)
    {
        std::string buffer;
        std::cin >> buffer;

        UInt number = 0;
        for (auto it = buffer.rbegin(); it != buffer.rend(); it++)
        {
            number *= 64;
            number += char_to_number(*it);
        }

        Zp x(number);
        message.push_back(EZp(x));
    }

    return message;
}

std::vector<std::pair<EZp, EZp>> encrypt_el_gamal(const EZp& s, const std::vector<EZp>& message)
{
    //std::mt19937 rnd(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    std::vector<std::pair<EZp, EZp>> result;
    for (const EZp &h_i : message)
    {
        //UInt b_i = UInt(rnd()) % (ord - 1) + 1;
        UInt b_i = 1;

        EZp first = g * b_i;
        assert(first == g);
        EZp second = h_i + (s * b_i);
        assert(second == (h_i + s));

        result.push_back({ first, second });
    }

    return result;
}

int main()
{
    /*std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.tie(0);*/

    /* Init global vars */

    EZp s;
    std::cin >> s;

    std::vector<EZp> message = scan_message();

    std::vector<std::pair<EZp, EZp>> result = encrypt_el_gamal(s, message);

    for (const auto &it : result)
    {
        std::cout << it.first << '\n';
        std::cout << it.second << '\n';
    }

    return 0;
}
