#include <bits/stdc++.h>

// https://github.com/dmkz/competitive-programming/blob/master/e-olymp.com/0317.cpp
const long double PI = std::acos(-1.0L);

struct UInt {
    static const int BASE = 1000000000; // Основание системы счисления
    static const int WIDTH = 9;       // Количество десятичных цифр, которые хранятся в одной цифре
    
    // Вектор под цифры числа:
    int digits[60];
    int sz;
    
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
    while (digits[sz - 1] == 0 && (int)sz > 1) sz--;
    return *this;
}   
 
// Конструктор от короткого целого
UInt::UInt(int64_t number) {
    sz = 0;
    do {
        digits[sz++] = number % BASE;
        number /= BASE;
    } while (number > 0);
    normalize();
}
 
// Конструктор от вектора из цифр:
UInt::UInt(const std::vector<int>& digits_) {
    sz = digits_.size();
    memcpy(digits, digits_.data(), sz * sizeof(int));
    normalize();
}
 
// Конструктор от строчки:
UInt::UInt(const std::string& s) {
    sz = 0;
    const int size = (int)s.size();
    for (int idGroup = 1, nGroups = size / WIDTH; idGroup <= nGroups; ++idGroup) {            
        digits[sz++] = std::stoi(s.substr(size-idGroup * WIDTH, WIDTH));
    }
    if (size % WIDTH != 0) {
        digits[sz++] = std::stoi(s.substr(0, size % WIDTH));
    }
    if (sz == 0) {
        digits[sz++] = 0;
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
        if (i >= (int)sz) digits[sz++] = 0;
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
    if (other.sz == 1u) {
        return *this += other.digits[0];
    }
    const int s1 = this->sz;
    const int s2 = other.sz;
    int rem = 0;
    for (int i = 0; i < s1 || i < s2 || rem > 0; ++i) {
        int d1 = i < s1 ? this->digits[i] : (digits[sz++] = 0, 0);
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
    for (int i = 0; i < (int)sz && rem < 0; ++i) {
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
    if (other.sz == 1u) {
        return *this -= other.digits[0];
    }
    const int s1 = this->sz;
    const int s2 = other.sz;
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
    for (int i = 0; i < sz; i++) {
        int& d = digits[i];
        rem += 1LL * d * num;
        auto div = rem / BASE;
        d = rem - div * BASE;
        rem = div;
    }
    if (rem > 0) digits[sz++] = rem;
    return this->normalize();
}
 
// Медленное произведение:
UInt UInt::mult(const UInt& other) const {
    if (other.sz == 1u) {
        return *this * other.digits[0];
    }
    const int s1 = (int)this->sz;
    const int s2 = (int)other.sz;
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
 
// Деление на короткое:
UInt& UInt::operator/=(const uint32_t num) {
    assert(num > 0);
    if (num >= BASE) {
        return *this /= UInt(num);
    }
    int64_t rem = 0;
    for (int j = (int)sz-1; j >= 0; --j) {
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
    for (int i = (int)a.sz-1; i >= 0; --i) {
        ((rem *= UInt::BASE) += a.digits[i]) %= num;
    }
    return rem;
}
 
// Целая часть и остаток от деления:
std::pair<UInt, UInt> UInt::div_mod(const UInt& other) const {
    if (other.sz == 1u) {
        return {std::move(*this / other.digits[0]), *this % other.digits[0]};
    }
    const int norm = BASE / (other.digits[other.sz - 1] + 1);
    const UInt a = *this * norm;
    const UInt b = other * norm;
    const int a_size = (int)a.sz;
    const int b_size = (int)b.sz;
    UInt q, r;
    if (q.sz > a_size)
    {
        q.sz = a_size;
    }
    while (q.sz < a_size)
    {
        q.digits[q.sz++] = 0;
    }
    for (int i = a_size - 1; i >= 0; --i) {
        r *= BASE;
        r += a.digits[i];
        int s1 = (int)r.sz <= b_size ? 0 : r.digits[b_size];
        int s2 = (int)r.sz <= b_size - 1 ? 0 : r.digits[b_size - 1];
        int d = (1LL * BASE * s1 + s2) / b.digits[b.sz - 1];
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
    if (this->sz > other.sz) return 1;
    if (this->sz < other.sz) return -1;
    for (int i = (int)sz-1; i >= 0; --i) {
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
    os << number.digits[number.sz - 1];
    for (int i = (int)number.sz-2; i >= 0; --i) {
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

UInt p("115792089210356248762697446949407573530086143415290314195533631308867097853951");
UInt p_minus_two("115792089210356248762697446949407573530086143415290314195533631308867097853949");
/**
 * Field Zp, p -- prime
 */
struct Zp
{
    UInt value;

    Zp(): value()
    {
    }
 
    Zp(const UInt& value_): value(value_)
    {
        /*if (value > p)
        {
            value %= p;
        }*/
    }
 
    Zp& operator*=(const Zp& other)
    {
        value *= other.value;
        if (value > p)
            value %= p;
        return (*this);
    }
 
    Zp& operator+=(const Zp& other)
    {
        value += other.value;
        if (value > p)
        {
            value -= p;
        }
        return (*this);
    }
 
    Zp& operator-=(const Zp& other)
    {
        if (value < other.value)
        {
            value += p;
        }
        value -= other.value;
        return (*this);
    }

    Zp& mul(const Zp& other)
    {
        value *= other.value;
        return (*this);
    }
 
    Zp& add(const Zp& other)
    {
        value += other.value;
        return (*this);
    }
 
    Zp& sub(const Zp& other)
    {
        value -= other.value;
        return (*this);
    }

    Zp& normalize()
    {
        value %= p;
        return (*this);
    }
};

Zp pow2(Zp a, int n)
{
    Zp res(1);
    while (n > 0)
    {
        if (n & 1)
        {
            res *= a;
        }
        a *= a;
        n >>= 1;
    }
    return res;
}

Zp pow(Zp a, const UInt& n)
{
    Zp res(1);
    for (size_t i = 0; i + 1 < n.sz; i++)
    {
        int x = n.digits[i];
        res *= pow2(a, x);
        a = pow2(a, UInt::BASE);
    }
    res *= pow2(a, n.digits[n.sz - 1]);
    return res;
}

Zp get_inv(const Zp& a)
{
    return pow(a, p_minus_two);
}

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
    Zp result(lhs);
    result += rhs;
    return result;
}
 
Zp operator-(const Zp& lhs, const Zp& rhs)
{
    Zp result(lhs);
    result -= rhs;
    return result;
}
 
Zp operator*(const Zp& lhs, const Zp& rhs)
{
    Zp result(lhs);
    result *= rhs;
    return result;
}

/* NIST P-256 */
Zp a(UInt("3")); // -3
Zp b(UInt("41058363725152142129326129780047268409114441015993725554835256314039467401291"));

struct EZp
{
    Zp x, y, z;
    bool is_inf;

    EZp(): x(), y(), z(), is_inf(false)
    {
    }

    EZp(const Zp& x_, const Zp& y_, const Zp& z_): x(x_), y(y_), z(z_), is_inf(false)
    {
        // not check
    }

    EZp& inverse()
    {
        if (!is_inf)
        {
            y = p - y;
        }
        return *this;
    }
};

EZp get_e()
{
    EZp tmp;
    tmp.is_inf = true;
    return tmp;
}

std::istream& operator>>(std::istream& is, EZp& ezp)
{
    std::string buffer;
    std::getline(is, buffer, '\n');
    if (buffer == "")
    {
        std::getline(is, buffer, '\n');
    }
    if ((buffer.size() == 1) && (buffer[0] == 'Z'))
    {
        ezp.is_inf = true;
    }
    else
    {
        std::string tmp;
        size_t i = 0;
        for (; i < buffer.size(); i++)
        {
            if (buffer[i] == ' ')
            {
                break;
            }
            tmp.push_back(buffer[i]);
        }
        ezp.x = Zp(UInt(tmp));
        tmp.resize(0);
        i++;
        for (; i < buffer.size(); i++)
        {
            tmp.push_back(buffer[i]);
        }
        ezp.y = Zp(UInt(tmp));
        ezp.z = Zp(1);
        ezp.is_inf = false;
    }
    return is;
}

Zp three(3);
Zp two(2);

EZp twice(const EZp& lhs)
{
    if ((lhs.is_inf) || (lhs.y == UInt(0)))
    {
        return get_e();
    }
    Zp tmp(a);
    tmp.mul(lhs.z);
    tmp.mul(lhs.z);
    tmp.normalize();
    Zp t = three;
    t.mul(lhs.x);
    t.mul(lhs.x);
    t.normalize();
    t -= tmp;
    Zp u = two;
    u.mul(lhs.y);
    u.mul(lhs.z);
    u.normalize();
    Zp v = two;
    v.mul(u);
    v.mul(lhs.x);
    v.mul(lhs.y);
    v.normalize();
    Zp w = t;
    w.mul(t);
    tmp = two;
    tmp *= v;
    w -= tmp;
    w.normalize();
    EZp ans(u, t, u);
    ans.x *= w;
    tmp = u;
    tmp.mul(lhs.y);
    tmp.mul(tmp);
    tmp.mul(two);
    tmp.normalize();
    ans.y.mul(v - w);
    ans.y.sub(tmp);
    ans.y.normalize();
    ans.z.mul(u);
    ans.z.mul(u);
    ans.z.normalize();
    return ans;
}

EZp operator+(const EZp& lhs, const EZp& rhs)
{
    const Zp& x0 = lhs.x;
    const Zp& y0 = lhs.y;
    const Zp& z0 = lhs.z;
    const Zp& x1 = rhs.x;
    const Zp& y1 = rhs.y;
    const Zp& z1 = rhs.z;
    if (lhs.is_inf)
    {
        return rhs;
    }
    if (rhs.is_inf)
    {
        return lhs;
    }
    Zp t0(y0);
    t0 *= z1;
    Zp t1(y1);
    t1 *= z0;
    Zp u0(x0);
    u0 *= z1;
    Zp u1(x1);
    u1 *= z0;
    if (u0 == u1)
    {
        if (t0 == t1)
        {
            return twice(lhs);
        }
        return get_e();
    }
    Zp t(t0);
    t -= t1;
    Zp u(u0);
    u -= u1;
    Zp u2(u);
    u2 *= u;
    Zp v(z0);
    v *= z1;
    Zp w(t);
    w.mul(t);
    w.mul(v);
    w.normalize();
    Zp tmp(u2);
    tmp *= u0 + u1;
    w -= tmp;
    Zp u3(u);
    u3 *= u2;
    EZp ans(u, t, u3);
    ans.x *= w;
    tmp = u0;
    tmp *= u2;
    tmp -= w;
    ans.y *= tmp;
    tmp = t0;
    tmp *= u3;
    ans.y -= tmp;
    ans.z *= v;
    return ans;
}

EZp pow2(EZp a, int n)
{
    EZp res = get_e();
    while (n > 0)
    {
        if (n & 1)
        {
            res = res + a;
        }
        a = twice(a);
        n >>= 1;
    }
    return res;
}

EZp operator*(EZp a, const UInt& n)
{
    EZp res = get_e();
    for (size_t i = 0; i + 1 < n.sz; i++)
    {
        int x = n.digits[i];
        res = res + pow2(a, x);
        a = pow2(a, UInt::BASE);
    }
    return res + pow2(a, n.digits[n.sz - 1]);
}

char number_to_char(int number)
{
    if (number == 62)
        return '_';
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

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.tie(0);

    UInt a = 0;
    std::cin >> a;

    int n = 0;
    std::cin >> n;

    for (int i = 0; i < n; i++)
    {
        EZp r_i, m_i;
        std::cin >> r_i >> m_i;

        m_i = m_i + (r_i * a).inverse();

        m_i.x *= get_inv(m_i.z);

        UInt number(m_i.x.value);

        while (number > 0)
        {
            std::cout << number_to_char(number % 64);
            number /= 64;
        }

        std::cout << '\n';
    }

    return 0;
}
