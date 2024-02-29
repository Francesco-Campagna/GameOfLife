#ifndef JALLEGRO_PANEL_H
#define JALLEGRO_PANEL_H

#include <allegro5/allegro.h>
#include "Graphics.h"
#include "Color.h"
#include <iostream>

namespace jallegro {
    class Panel {
    public:
        Panel() {
            backgroundColor = Color::white();
            closed_ = false;
        }

        virtual ~Panel() {}

        void repaint() {
            if(closed_) return;
            this->paintComponent(g);
            al_flip_display();
        }

        void setBackgroundColor(Color c) {
            backgroundColor = c;
        }

        void close() {
            closed_ = true;
        }

        bool closed() {
            return closed_;
        }

    protected:
        virtual void paintComponent(Graphics g) {
            al_clear_to_color(backgroundColor.toAllegro());
        }

    private:
        Graphics g;
        Color backgroundColor;
        bool closed_;
    };
}

#endif