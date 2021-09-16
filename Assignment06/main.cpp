#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <vector>

using std::vector;
using std::logic_error;
void LogicError(const char* format, ...) {
    char buffer[4096];
    va_list args;
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);
    throw logic_error(buffer);
}

class Vertex {
  private:
    double x_;
    double y_;
    size_t index_;

  public:
    Vertex(const double x, const double y, const size_t index)
    : x_(x), y_(y), index_(index) { }
    double x() const {
        return x_;
    }
    double y() const {
        return y_;
    }
    size_t index() const {
        return index_;
    }
    static double distance(const Vertex& a, const Vertex& b) {
        double dx = a.x_ - b.x_;
        double dy = a.y_ - b.y_;
        return sqrt(dx * dx + dy * dy);
    }
} ;

class Chord {
  private:
    const Vertex* one_;
    const Vertex* two_;
    double dist_;

  public:
    Chord() : one_(NULL), two_(NULL), dist_(INFINITY) { }
    Chord(const Vertex& one, const Vertex& two)
    : one_((one.index() < two.index()) ? (&one) : (&two)),
      two_((one.index() < two.index()) ? (&two) : (&one)),
      dist_(Vertex::distance(one, two)) { }
    double dist() const {
        return dist_;
    }
    const Vertex* one() const {
        return one_;
    }
    const Vertex* two() const {
        return two_;
    }
    void print() const {
        printf("%3lu:(%12.4f,%12.4f) => %3lu:(%12.4f,%12.4f) dist: %12.4f\n",
               one_->index(), one_->x(), one_->y(),
               two_->index(), two_->x(), two_->y(),
               dist_);
    }
    bool operator ==(const Chord& other) const {
        return (one_->index() == other.one()->index() &&
                two_->index() == other.two()->index());
    }
} ;

class Polygon {
  private:
    const char*     path_;
    vector<Vertex>  points_;
    vector<Chord>   chords_;
    double**        value_matrix_;
    double**        dists_matrix_;
    vector<Chord>** chord_matrix_;
    void read() {
        FILE* file;
        if ((file = fopen(path_, "r")) == NULL)
            LogicError("fopen('%s', 'r') error: %s", path_, strerror(errno));
        double x, y;
        while (fscanf(file, "%lf%lf", &x, &y) == 2)
            points_.push_back(Vertex(x, y, points_.size()));
        fclose(file);
    }
    void add_chord(const Chord& c, vector<Chord> * const out) {
        // Don't add if edge.
        if ((c.one()->index() + 1) % points_.size() == c.two()->index() ||
            (c.two()->index() + 1) % points_.size() == c.one()->index())
            return;
        // Don't add if chord exists in chords_
        for (auto i = out->begin(); i != out->end(); ++i)
            if (c == *i) return;
        out->push_back(c);
    }
    static void copy_to(vector<Chord> * const src, vector<Chord> * const dst) {
        for (auto i = src->begin(); i != src->end(); ++i)
            dst->push_back(*i);
    }
    double triangle_perim(size_t i, size_t j, size_t k) {
        return dists_matrix_[i][j] + dists_matrix_[j][k] +
               dists_matrix_[k][i];
    }
    double decompose(const size_t i, const size_t j, vector<Chord>* const out) {
        if (value_matrix_[i][j] >= 0) {
            copy_to(&chord_matrix_[i][j], out);
            return value_matrix_[i][j];
        } else if (j == i + 1) {
            return value_matrix_[i][j] = 0;
        } else {
            double min = INFINITY;
            size_t best_k;
            vector<Chord> best_out;
            for (size_t k = i + 1; k < j; ++k) {
                vector<Chord> new_out;
                double v = decompose(i, k, &new_out)
                         + decompose(k, j, &new_out)
                         + triangle_perim(i, j, k);
                if (v < min) {
                    min = v;
                    best_k = k;
                    best_out = new_out;
                }
            }
            add_chord(Chord(points_[i],      points_[j]),      &best_out);
            add_chord(Chord(points_[j],      points_[best_k]), &best_out);
            add_chord(Chord(points_[best_k], points_[i]),      &best_out);
            copy_to(&best_out, out);
            chord_matrix_[i][j]   = best_out;
            return value_matrix_[i][j] = min;
        }
    }

  public:
    explicit Polygon(const char* path) : path_(path) {
        read();
        value_matrix_ = new double*[points_.size()];
        dists_matrix_ = new double*[points_.size()];
        chord_matrix_ = new vector<Chord>*[points_.size()];
        for (size_t i = 0; i < points_.size(); ++i) {
            value_matrix_[i] = new double[points_.size()];
            dists_matrix_[i] = new double[points_.size()];
            chord_matrix_[i] = new vector<Chord>[points_.size()];
            for (size_t j = 0; j < points_.size(); ++j) {
                value_matrix_[i][j] = -1;
                dists_matrix_[i][j] = Vertex::distance(points_[i], points_[j]);
            }
        }
    }
    double decompose() {
        return decompose(0, points_.size() - 1, &chords_);
    }
    void print_chords() {
        for (auto i = chords_.begin(); i != chords_.end(); ++i)
            i->print();
    }
} ;
int main() {
    Polygon p("polygon2.txt");
    printf("size: %f\n", p.decompose());
    p.print_chords();
}
