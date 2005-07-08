/** "C" event display.
 * Main part. 
 *
 * Alexey Zhelezov, DESY/ITEP, 2005 */

/**
 * This file is internal. It must not be
 * included into enduser application.
 */


#ifndef __CED_H
#define __CED_H

typedef void (*ced_draw_cb)(void *data);


/**
 * Register new element type. Order is important!
 * Appropriate code must be defined in both
 * client and server parts.
 *
 *  item_size - size of one item in bytes.
 *  draw_func - function to call to draw one item,
 *              called from ced_do_draw_event()
 *              not used on client side.
 */
unsigned ced_register_element(unsigned item_size,ced_draw_cb draw_func);

/**
 * To be called from element functions
 * on client side.
 * Return allocated space for one item
 * with size, specified by ced_register_element()
 *
 * Example: assume struct Dummy { int i; }; is item.
 *
 *          DummyID=ced_register_element(sizeof(struct Dummy),0);
 *          ...
 *          struct Dummy *item=(struct Dummy *)ced_add_element(DummyID);
 *            item->i=0;
 * This will add one Dummy item to the event
 */
void *ced_add(unsigned id);

/**
 * To be called in paint function
 *
 * It calls user defined functions for
 * each item of all elements types.
 */
void ced_do_draw_event(void);

/**
 * Server side function.
 * Must be used to process all incoming
 * messages from client.
 *
 * It return positive value when
 * new event must be drawn.
 *
 * Example:
 *      glut_tcp_server(7285,my_process_input)
 *
 *      my_process_input(x){
 *        if(ced_process_input(x)>0)
 *          <do redraw>
 */
int ced_process_input(void *data);

#endif /* __CED_H  */
