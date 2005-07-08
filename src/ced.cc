/* "C" event display.
 * Communications related part. 
 *
 * Alexey Zhelezov, DESY/ITEP, 2005 
 *
 * F.Gaede, DESY : made compliant with g++ ( introducted casts where necessary ) 
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>

#include <ced.h>

static int ced_fd=-1; // CED connection socket

static unsigned short ced_port=7927; // port No of CED (assume localhost)

// Return 0 if can be connected, -1 otherwise.
static int ced_connect(void){
  static time_t last_attempt=0;
  time_t ct;
  struct sockaddr_in addr;

  if(ced_fd>=0)
    return 0; // already connected;
  time(&ct);
  if(ct-last_attempt<5)
    return -1; // don't try reconnect all the time
  addr.sin_family=AF_INET;
  addr.sin_port=htons(ced_port);
  addr.sin_addr.s_addr=htonl(0x7f000001); // 127.0.0.1
  ced_fd=socket(PF_INET,SOCK_STREAM,0);
  if(connect(ced_fd,(struct sockaddr *)&addr,sizeof(addr))){
    if(!last_attempt)
      perror("WARNING:CED: can't connect to CED");
    time(&last_attempt);
    close(ced_fd);
    ced_fd=-1;
    return -1;
  }
  fprintf(stderr,"INFO:CED: connected to CED\n");
  return 0;
}


typedef struct {
  unsigned size;     // size of one item in bytes
  unsigned char *b;  // "body" - data are stored here
                     // (here is some trick :)
  unsigned count;    // number of usefull items
  unsigned alloced;  // number of allocated items
  ced_draw_cb draw;  // draw fucation, NOT used in CED client
} ced_element;

typedef struct {
  ced_element *e;
  unsigned      e_count;
} ced_event;

static ced_event eve = {0,0};
// NOT used in CED client
static ced_event ceve = {0,0}; // current event on screen

// we reserve this size just before ced_element.b data
#define HDR_SIZE 8

unsigned ced_register_element(unsigned item_size,ced_draw_cb draw_func){
  ced_element *pe;
  if(!(eve.e_count&0xf))
    eve.e=(ced_element* ) realloc(eve.e,(eve.e_count+0x10)*sizeof(ced_element));
  pe=eve.e+eve.e_count;
  memset(pe,0,sizeof(*pe));
  pe->size=item_size;
  pe->draw=draw_func;
  return eve.e_count++;
}

static void ced_reset(void){
  unsigned i;
  for(i=0;i<eve.e_count;i++)
    eve.e[i].count=0;
}

static void ced_buf_alloc(ced_element *pe,unsigned count){
  if(!pe->b)
    pe->b=(unsigned char*) malloc(count*pe->size+HDR_SIZE);
  else
    pe->b=(unsigned char*) realloc(pe->b-HDR_SIZE,count*pe->size+HDR_SIZE);
  pe->b+=HDR_SIZE;
  pe->alloced=count;
}

void *ced_add(unsigned id){
  ced_element *pe;
  if(id >= eve.e_count){
    fprintf(stderr,"BUG:CED: attempt to access not registered element\n");
    return 0;
  }
  pe=eve.e+id;
  if(pe->count==pe->alloced)
    ced_buf_alloc(pe,pe->alloced+256);
  return (pe->b+(pe->count++)*pe->size);
}

static void ced_event_copy(ced_event *trg){
  unsigned i;
  ced_element *pe;
  if(trg->e_count<eve.e_count)
    trg->e= (ced_element*) realloc(trg->e,eve.e_count*sizeof(ced_element));
  for(i=0;i<eve.e_count;i++){
    pe=trg->e+i;
    if(i<trg->e_count){
      if(pe->alloced<eve.e[i].alloced)
	ced_buf_alloc(pe,eve.e[i].alloced);
      pe->count=eve.e[i].count;
    } else {
      memcpy(pe,eve.e+i,sizeof(ced_element));
      if(pe->b){
	pe->b=0;
	ced_buf_alloc(pe,pe->alloced);
      }
    }
    if(pe->count)
      memcpy(pe->b,eve.e[i].b,pe->count*pe->size);
  }
  trg->e_count=eve.e_count;
}

void ced_do_draw_event(void){
  unsigned i,j;
  ced_element *pe;
  unsigned char *pdata;
  for(i=0;i<ceve.e_count;i++){
    pe=ceve.e+i;
    if(!pe->draw)
      continue;
    for(pdata=pe->b,j=0;j<pe->count;j++,pdata+=pe->size)
      (*(pe->draw))(pdata);
  }
}

typedef enum {
  DRAW_EVENT=10000
} MSG_TYPE;

int ced_process_input(void *data){
  struct _phdr{
    unsigned size;
    unsigned type;
    unsigned char b[4];
  } *hdr = (_phdr*) data;
  unsigned count;
  ced_element *pe;
  
  if(!data){ // new client is connected
    ced_reset();
    return 0;
  }

  if(hdr->type == DRAW_EVENT){
    ced_event_copy(&ceve);
    ced_reset();
    return 1;
  }
  if(hdr->type>=eve.e_count){
    fprintf(stderr,"WARNING:CED: undefined element type (%u), ignored\n",
	    hdr->type);
    return 0;
  }
  pe=eve.e+hdr->type;
  if((hdr->size-HDR_SIZE)%pe->size){
    fprintf(stderr,"BUG:CED: size alignment is wrong for element %u\n",
	    hdr->type);
    return 0;
  }
  count=(hdr->size-HDR_SIZE)/pe->size;
  if(!count)
    return 0;
  if(count>=pe->alloced)
    ced_buf_alloc(pe,count+256);
  memcpy(pe->b,hdr->b,count*pe->size);
  pe->count=count;
  return 0;
}

void ced_send_event(void){
  struct _phdr{
    unsigned size;
    unsigned type;
  } *hdr,draw_hdr;
  unsigned i,sent_sum,problem=0;
  char *buf;
  int sent;
  ced_element *pe;

  if(ced_connect())
    return;
  for(i=0;i<eve.e_count && !problem;i++){
    pe=eve.e+i;
    if(!pe->count)
      continue;
    hdr=(struct _phdr *)(pe->b-HDR_SIZE); // !!! HERE is the trick :)
    hdr->type=i;
    hdr->size=HDR_SIZE+pe->count*pe->size;
    sent_sum=0;
    buf=(char *)hdr;
    while(sent_sum<hdr->size){
	sent=write(ced_fd,buf+sent_sum,hdr->size-sent_sum);
	if(sent<0){
	    problem=1;
	    break;
	}
	sent_sum+=sent;
    }
  }
  if(!problem){
    draw_hdr.size=HDR_SIZE;
    draw_hdr.type=DRAW_EVENT;
    if(write(ced_fd,&draw_hdr,HDR_SIZE)!=HDR_SIZE)
      problem=1;
  }
  if(problem){
    perror("WARNING:CED: can't send event, till next time...");
    close(ced_fd);
    ced_fd=-1;
  }
}

#include <signal.h>
// API
void ced_client_init(const char *host,unsigned short port){
  ced_port=port;
  signal(SIGPIPE,SIG_IGN);
}

void ced_new_event(void){
  ced_reset();
}

void ced_draw_event(void){
  ced_send_event();
  ced_reset();
}
