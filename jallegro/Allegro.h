#ifndef JALLEGRO_ALLEGRO_H
#define JALLEGRO_ALLEGRO_H

#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <iostream>
namespace jallegro {
    class Allegro {
    public:
        static void init() {
            if (!al_init() || !al_init_primitives_addon()) {
                std::cerr << "Unable to init allegro!" << std::endl;
                exit(10);
            }
        }

        static void sleep(float seconds) {
            al_rest(seconds);
        }
    };
}

#endif