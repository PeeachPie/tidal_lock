#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include <stdexcept>
#include "constants.hpp"
#include "particle.hpp"

struct Data {
    glm::dvec2 sum_xm_{ 0.0, 0.0 }; // sum x_i * m_i
    double total_mass_ = 0.0;
    int id_ = -1;

    void clear_id() {
        id_ = -1;
    }

    void append_particle(Particle p, double p_mass) {
        sum_xm_ += p.pos * p_mass;
        total_mass_ += p_mass;

        clear_id();
    }

    void add_particle(Particle p, double p_mass, int id) {
        append_particle(p, p_mass);
        id_ = id;
    }

    glm::dvec2 mass_center() {
        if (total_mass_ == 0.0)
            throw std::runtime_error("`mass_center` called on empty Node");
        return sum_xm_ / total_mass_;
    }
};


class Node {
public: // TODO: private
    // t - top, b - bottom
    Node *tl_ = nullptr; 
    Node *tr_ = nullptr; 
    Node *bl_ = nullptr;
    Node *br_ = nullptr;

    glm::dvec2 tl_bound_, br_bound_;

private:
    bool empty_;
    glm::dvec2 point_pos_;
    Data data_;

public:
    Node(glm::dvec2 tl_bound, glm::dvec2 br_bound): 
        tl_bound_(tl_bound), br_bound_(br_bound), empty_(true) {};

    Node(): tl_bound_({-INF, INF}), br_bound_({INF, -INF}), empty_(true) {};

    void split() {
        assert(std::min(std::fabs(br_bound_.x - tl_bound_.x), std::fabs(br_bound_.y - tl_bound_.y)) > 10);

        if (is_splitted())
            throw std::runtime_error("`split` on splitted Node");

        glm::dvec2 mid = (tl_bound_ + br_bound_) / 2.0;

        tl_ = new Node(tl_bound_, mid);
        tr_ = new Node({ {mid.x, tl_bound_.y}, {br_bound_.x, mid.y} });
        bl_ = new Node({ {tl_bound_.x, mid.y}, {mid.x, br_bound_.y} });
        br_ = new Node(mid, br_bound_);

        Node *new_place;

        if (tl_->in_bounds(point_pos_)) new_place = tl_;
        else if (tr_->in_bounds(point_pos_)) new_place = tr_;
        else if (bl_->in_bounds(point_pos_)) new_place = bl_;
        else if (br_->in_bounds(point_pos_)) new_place = br_;

        push(new_place, point_pos_, data_);
    }

    void push(Node *new_place, glm::dvec2 pos, Data data) {
        new_place->move_particle(pos, data);
        data_.clear_id();
    }

    void move_particle(glm::dvec2 pos, Data data) {
        point_pos_ = pos;
        data_ = data;
        empty_ = false;
    }

    void append_particle(Particle p, double p_mass) {
        data_.append_particle(p, p_mass);
    }

    void add_particle(Particle p, double p_mass, int id) {
        point_pos_ = p.pos;
        data_.add_particle(p, p_mass, id);
        empty_ = false;
    }

    glm::dvec2 get_point_pos() const { // TODO: if empty
        return point_pos_;
    }

    Data get_data() const { // TODO: if empty
        return data_;
    }

    bool empty() const {
        return empty_;
    }

    bool is_splitted() const {
        return tl_ != nullptr && tr_ != nullptr && br_ != nullptr && bl_ != nullptr; 
    }

    bool in_bounds(glm::dvec2 pos) const {
        return pos.x <= br_bound_.x && pos.x > tl_bound_.x && 
               pos.y >= br_bound_.y && pos.y < tl_bound_.y;
    }

    ~Node() {
        delete tl_;
        delete tr_;
        delete bl_;
        delete br_;
    }
};

class BarnesHutQuadTree {
private:
    Node *root_ = nullptr;
public:
    int cnt_ = 0; // TODO: remove
    double theta_;

    BarnesHutQuadTree(double theta=0.5): root_(new Node()), theta_(theta) {};

    BarnesHutQuadTree(double bound, double theta): 
        root_(new Node({-bound, bound}, {bound, -bound})),
        theta_(theta) 
    {};

    void insert(Particle &p, double p_mass, int id) {
        assert(root_->in_bounds(p.pos));
        _insert(root_, p, p_mass, id);
    }

    // TODO: id only
    glm::dvec2 calc_force(Particle &p, double p_mass, int id) const {
        assert(root_->in_bounds(p.pos));
        return _calc_force(root_, p, p_mass, id);
    }

    void print() const {
        print_node(root_, 0);
    }

    ~BarnesHutQuadTree() {
        delete root_;
    }
    
private:
    void _insert(Node *node, Particle &p, double p_mass, int id) {
        if (node == nullptr)
            throw std::runtime_error("`_insert` called on null");

        if (!node->in_bounds(p.pos)) return;

        if (node->empty()) {
            node->add_particle(p, p_mass, id);
        }
        else {
            if (!node->is_splitted()) {
                node->split();
                cnt_ += 4;
            }
            node->append_particle(p, p_mass);

            // TODO: -> if
            _insert(node->tl_, p, p_mass, id);
            _insert(node->tr_, p, p_mass, id);
            _insert(node->br_, p, p_mass, id);
            _insert(node->bl_, p, p_mass, id);
        }
    }

    glm::dvec2 _calc_force(Node *node, Particle &p, double p_mass, int id) const {
        Data data = node->get_data();

        if (data.id_ == id || node->empty()) {
            return { 0.0, 0.0 };
        }

        glm::dvec2 mass_center = data.mass_center();

        double d = glm::distance(mass_center, p.pos);
        double s = node->br_bound_.x - node->tl_bound_.x;

        if (s/d > theta_ && node->is_splitted()) {
            return 
                _calc_force(node->tl_, p, p_mass, id) + 
                _calc_force(node->tr_, p, p_mass, id) +
                _calc_force(node->bl_, p, p_mass, id) + 
                _calc_force(node->br_, p, p_mass, id);
        }
        else {
            // std::cout << "debug: " << data.id_ << ' ' << data.total_mass_ << '\n';

            double r = std::max(glm::distance(p.pos, mass_center), SMOOTHING_LENGTH * 0.3); // TODO: explicit
            glm::dvec2 n = (mass_center - p.pos) / r;
            glm::dvec2 f = (G * (data.total_mass_ * p_mass) / (r * r)) * n;

            return f;
        }
    }

    void print_node(const Node* node, int depth) const {
        if (node == nullptr)
            return;

        std::string indent(depth * 4, ' ');

        std::cout << indent << "Node ";

        std::cout << "[(" << node->tl_bound_.x << ", " << node->tl_bound_.y
                << ") -> ("
                << node->br_bound_.x << ", " << node->br_bound_.y << ")] ";

        if (node->empty()) {
            std::cout << "EMPTY";
        } else {
            const Data data = node->get_data();

            std::cout
                << "mass=" << data.total_mass_
                << "  sum_xm=(" << data.sum_xm_.x << ", " << data.sum_xm_.y << ")";

            if (data.id_ != -1)
                std::cout << "  particle_id=" << data.id_;
            else
                std::cout << "  particle_id=<internal>";

            auto p = node->get_point_pos();
            std::cout << "  pos=(" << p.x << ", " << p.y << ")";
        }

        if (node->is_splitted())
            std::cout << " SPLIT";

        std::cout << '\n';

        if (node->is_splitted()) {
            std::cout << indent << " TL:\n";
            print_node(node->tl_, depth + 1);

            std::cout << indent << " TR:\n";
            print_node(node->tr_, depth + 1);

            std::cout << indent << " BL:\n";
            print_node(node->bl_, depth + 1);

            std::cout << indent << " BR:\n";
            print_node(node->br_, depth + 1);
        }
    }
};

