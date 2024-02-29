#ifndef JALLEGRO_FRAME_H
#define JALLEGRO_FRAME_H

#include <allegro5/allegro.h>
#include "Panel.h"
#include <iostream>
using namespace std;

namespace jallegro {

    void* Func_Thread(ALLEGRO_THREAD* thr, void* arg);

    class Frame {
    public:
        Frame(unsigned int width, unsigned int height, const char* title) {
            display = al_create_display(width, height);
            if (display == nullptr)
                return;
            al_set_window_title(display, title);
            al_clear_to_color(Color::white().toAllegro());
            setExitOnClose(true);
            shown = false;
        }

        virtual ~Frame() {
            if(display != nullptr)
                al_destroy_display(display);
        }

        void show() {
            if(shown)
                return;
            shown = true;
            ALLEGRO_THREAD *thread_1 = al_create_thread(Func_Thread, this);
            al_start_thread(thread_1);
            if(panel != nullptr)
                panel->repaint();
        }

        void add(Panel *panel) {
            this->panel = panel;
        }

        void close() {
            if(panel != nullptr)
                panel->close();
        }

        void setExitOnClose(bool value) {
            exit = value;
        }

        bool exitOnClose() {
            return exit;
        }

    private:
        ALLEGRO_DISPLAY* display;
        Panel *panel;
        bool exit;
        friend void* Func_Thread(ALLEGRO_THREAD* thr, void* arg);
        bool shown;
    };

    void* Func_Thread(ALLEGRO_THREAD* thr, void* arg) {
        ALLEGRO_EVENT_QUEUE* event_queue = al_create_event_queue();
        if (!event_queue)
            exit(-1);
        Frame* f = (Frame*) arg;
        al_register_event_source(event_queue, al_get_display_event_source(f->display));
        ALLEGRO_EVENT ev;
        while (true) {
            al_wait_for_event(event_queue, &ev);
            if (ev.type == ALLEGRO_EVENT_DISPLAY_CLOSE)
                break;
        }
        al_destroy_display(f->display);
        f->display = nullptr;
        if(f->exitOnClose())
            exit(0);
        f->close();
        return NULL;
    }
}

#endif