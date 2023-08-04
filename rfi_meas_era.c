/* RFI_MEAS - Acquire trace from spectrum analyzer through LAN (TCP) */
/* Santosh */
/*
 RFI_MEAS creates a socket, establishes TCP connection using port 5025
and given static IP address. It then creates a header file to store the
values of all important parameters such as start and stop frequency, sweep
time and points, etc. Later it initiates a single trigger, collects the
trace and stores it in a .CSV file. The number of traces to sweep has to
be given as an argument while executing the program.
*/
/*
c***************************************************************
c History
c
c Santosh 12Nov14 original version
c Santosh 17Nov14 modified for multiple trace retrieval
c Santosh 27Jan15 modified to apply user settings and write the trace as CSV
c***************************************************************
*/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <bits/socket.h>
#include <string.h>
#include <sys/select.h>
#include <sys/time.h>
#include <sys/timeb.h>
#include <errno.h>
#include <netinet/tcp.h>
#include <netinet/in.h>
#include <time.h>

#define SOCKETS_IP_ADDRESS			"192.168.1.97" 
#define SOCKETS_PORT				5025
#define SOCKETS_BUFFER_SIZE			1024
#define SOCKETS_TIMEOUT				20

const struct tm *tm ;
struct tm *cur_time_struct;
struct timeb t_current;
unsigned short millitime;
float timesec;
long cur_time_long;
int time_info1[7];
char outfile[70],outhdr[70];
	

int WriteString(int MySocket,char string[])
{
	int retval;
	
	if((retval=send(MySocket,string,strlen(string),0))==-1) {
		printf("Error: Unable to send message (%i)...\n",errno);
		perror("sockets"); // Print error message based on errno
		exit(1);
	}
	return retval;
}

int WaitForData(int MySocket)
{
	fd_set MyFDSet;
	struct timeval tv;
	int retval;
	
	// Wait for data to become available
	FD_ZERO(&MyFDSet); 					// Initialize fd_set structure
	FD_SET(MySocket,&MyFDSet); 				// Add socket to "watch list"
	tv.tv_sec=SOCKETS_TIMEOUT; tv.tv_usec=0; 		// Set timeout
	retval=select(MySocket+1,&MyFDSet,NULL,NULL,&tv); 	// Wait for change
	
	// Interpret return value
	if(retval==-1) {
		printf("Error: Problem with select (%i)...\n",errno);
		perror("sockets"); // Print error message based on errno
		exit(1);
	}
	return retval; // 0 = timeout, 1 = socket status has changed
}

int ReadString(int MySocket,char *buffer)
{
	size_t actual;
	
	// Wait for data to become available
	if(WaitForData(MySocket)==0)
	{
		// Timeout
		printf("Warning: Query Timeout \n");
		return 0;
	}
	
	// Read data
	if((actual=recv(MySocket,buffer,100,0))==-1) {
		printf("Error: Unable to receive data (%i)...\n",errno);
		perror("sockets"); // Print error message based on errno
		exit(1);
	}
	else {
	buffer[actual]=0;
	}
	return actual;
}

void Setup(int MySocket)
{
	char SocketsBuffer[SOCKETS_BUFFER_SIZE];
	
	// Clear status and reset instrument to user settings
	WriteString(MySocket,"INST 'SA';*OPC?\n");
	ReadString(MySocket,SocketsBuffer);	
	WriteString(MySocket,"FREQ:STAR 500E6;*WAI;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);
	WriteString(MySocket,"FREQ:STOP 6.5E9;*WAI;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);
	WriteString(MySocket,":SENS:SWE:POIN 1001;:SENS:SWE:TYPE AUTO;*WAI;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);
	WriteString(MySocket,":SENS:BAND:AUTO 0;:SENS:BAND 600E3;*WAI;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);
	WriteString(MySocket,":SENS:BAND:VID:AUTO 0;:SENS:BAND:VID 1.5E6;*WAI;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);
	WriteString(MySocket,":SENS:POW:GAIN 0;:SENS:POW:ATT:AUTO 0;:SENS:POW:ATT 0;*WAI;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);
	WriteString(MySocket,"DISP:WIND:TRAC1:Y:RLEV 0;*WAI;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);
	WriteString(MySocket,"AMPL:SCAL LOG;*WAI;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);

	printf("User settings set\n");
}

void ReadTrace(int MySocket,char *outfile,int adc)
{
	int i,j,lcount,trace_size=0,buf_size=0;
	char *ctoken;
	char buffer[2000],trace_data[17000];
	int sec,msec;
	struct timeb current;
	FILE *fp;
	size_t length;
	
	time(&cur_time_long);
  	cur_time_struct = (struct tm *) localtime(&cur_time_long);
	ftime(&t_current);
	time_info1[0]=cur_time_struct->tm_hour;
	time_info1[1]=cur_time_struct->tm_min;
	time_info1[2]=cur_time_struct->tm_sec;
	time_info1[3]=cur_time_struct->tm_mday;
	time_info1[4]=cur_time_struct->tm_mon+1;
	time_info1[5]=cur_time_struct->tm_year + 1900;
	millitime = t_current.millitm;
	timesec = (float)cur_time_struct->tm_sec + (float)millitime/1000.0;
	
	// Wait for data to become available
	if(WaitForData(MySocket)==0)
	{
		// Timeout
		printf("Warning: Trace Timeout \n");
		Setup(MySocket);
		return;
	}
	
	ftime(&current);
	while(trace_size < 16010) // Change the iteration condition according to the sweep points
    {
        memset(buffer, 0,sizeof(buffer));  //clear the variable
		
		// Read data
		/*struct timeval tv;
		tv.tv_sec = 10;  // 10 Secs Timeout
		if(setsockopt(MySocket, SOL_SOCKET, SO_RCVTIMEO,(struct timeval *)&tv,sizeof(struct timeval)))
		{
			printf("\n Timeout!");
			break;
		}*/
		buf_size=recv(MySocket,buffer,sizeof(buffer),0);
		memcpy(trace_data+trace_size,buffer,buf_size);	
		trace_size += buf_size;
		//printf("%d %d\n",buf_size,trace_size);
	}
	trace_data[trace_size] = 0;
	//printf("%s",trace_data);
	if (trace_size >= 16010)
	{
		sec = (current.time % 86400) + 19800;
		msec = current.millitm;
		fp = fopen(outfile,"a");
		fprintf(fp,"%d,%d.%.01d,",adc,sec,msec);
		fprintf(fp,"%s",trace_data);
		fclose(fp); 
	}
	trace_size = 0;
	return;
}

void DeviceClear(int MySocket,char *buffer)
{
	WriteString(MySocket,"DCL\n");
	if(ReadString(MySocket,buffer)==0)
		return;
	if(strcmp(buffer,"DCL\n")==0)
		printf("DCL\\n received back from instrument...\n");
	else {
		printf("Warning: DCL response: %s\n",buffer);
	}
	return;
}

void SetNODELAY(int MySocket)
{
	int StateNODELAY = 1; // Turn NODELAY on
	int ret;
	
	ret=setsockopt(MySocket, 	// Handle to socket connection
		IPPROTO_TCP, 		// Protocol level (TCP)
		TCP_NODELAY, 		// Option on this level (NODELAY)
		(void *)&StateNODELAY, 	// Pointer to option variable
		sizeof StateNODELAY); 	// Size of option variable
	
	if(ret==-1) {
		printf("Error: Unable to set NODELAY option (%i)...\n",errno);
		perror("sockets"); // Print error message based on errno
		exit(1);
	}
	return;
}

int main(int argc,char **argv)
{
	int MySocket,MyControlPort;
	struct hostent *hp;
	char SocketsBuffer[SOCKETS_BUFFER_SIZE];
	char *trace_buf;
	static char *ctoken;
	struct sockaddr_in MyAddress,MyControlAddress;
	unsigned int ControlPort;
	int status,lcount,i,retval,tcount=0,adc;
	FILE *fh;
	clock_t begin, end;
	double time_spent;
	struct timeval timeout;      
    timeout.tv_sec = 10;
    timeout.tv_usec = 0;


	time(&cur_time_long);
  	cur_time_struct = (struct tm *) localtime(&cur_time_long);
	time_info1[0]=cur_time_struct->tm_hour;
	time_info1[1]=cur_time_struct->tm_min;
	time_info1[2]=cur_time_struct->tm_sec;
	time_info1[3]=cur_time_struct->tm_mday;
	time_info1[4]=cur_time_struct->tm_mon+1;
	time_info1[5]=cur_time_struct->tm_year + 1900;
	sprintf(outfile,"%04d-%02d-%02d_%02d%02d%02d.csv%c",
		  time_info1[5],time_info1[4],time_info1[3],
		  time_info1[0],time_info1[1],time_info1[2],'\0');
	sprintf(outhdr,"%04d-%02d-%02d_%02d%02d%02d_hdr%c",
		  time_info1[5],time_info1[4],time_info1[3],
		  time_info1[0],time_info1[1],time_info1[2],'\0');

	// Create socket (allocate resources)
	if((MySocket=socket(
		PF_INET, // IPv4
		SOCK_STREAM, // TCP
		0))==-1) {
 		printf("Error: Unable to create socket (%i)...\n",errno);
		perror("sockets"); // Print error message based on errno
		exit(1);
	}

	// Establish TCP connection
	memset(&MyAddress,0,sizeof(struct sockaddr_in)); 		// Set structure to zero
	MyAddress.sin_family=PF_INET; 					// IPv4
	MyAddress.sin_port = htons(SOCKETS_PORT); 			// Port number (in network order)
	MyAddress.sin_addr.s_addr = inet_addr(SOCKETS_IP_ADDRESS); 	// IP address (in network order)
	setsockopt(MySocket, SOL_SOCKET, SO_RCVTIMEO, (char *)&timeout, sizeof(timeout));
	
	if(connect(MySocket,(struct sockaddr *)&MyAddress,sizeof(struct sockaddr_in))==-1) 
	{
		printf("Error: Unable to establish connection to socket (%i)...\n",errno);
		perror("sockets"); // Print error message based on errno
		exit(1);
	}
	// Minimize latency by setting TCP_NODELAY option
	SetNODELAY(MySocket);
	Setup(MySocket);

	fh=fopen(outhdr,"w");
	
	// Get instrument's ID string
	WriteString(MySocket,"*IDN?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);
	fprintf(fh,"Instrument ID: %s\n",SocketsBuffer);
	
	// Get system date
	WriteString(MySocket,"SYST:DATE?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);
	fprintf(fh,"System date: %s",SocketsBuffer);
	
	// Get start and stop frequency
	WriteString(MySocket,":SENS:FREQ:STAR?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);
	fprintf(fh,"Start frequency: %s",SocketsBuffer);
	WriteString(MySocket,":SENS:FREQ:STOP?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);
	fprintf(fh,"Stop frequency: %s",SocketsBuffer);

	// Get center frequency
	WriteString(MySocket,":SENS:FREQ:CENT?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);
	fprintf(fh,"Center frequency: %s",SocketsBuffer);
	
	// Get frequency span
	WriteString(MySocket,":SENS:FREQ:SPAN?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);
	fprintf(fh,"Frequency span: %s",SocketsBuffer);
	
	// Get Res BW and Video BW values
	WriteString(MySocket,":SENS:BAND:RES?\n");
	ReadString(MySocket,SocketsBuffer);
	fprintf(fh,"Resolution BW: %s",SocketsBuffer);
	WriteString(MySocket,":SENS:BAND:VID?\n");
	ReadString(MySocket,SocketsBuffer);
	fprintf(fh,"Video BW: %s",SocketsBuffer);
	
	// Get the no. of sweep points
	WriteString(MySocket,":SENS:SWE:POIN?\n");
	ReadString(MySocket,SocketsBuffer);
	fprintf(fh,"Sweep points: %s",SocketsBuffer);
	
	// Get sweep time
	WriteString(MySocket,":SENS:SWE:ACQ?\n");
	ReadString(MySocket,SocketsBuffer);
	fprintf(fh,"Sweep time: %s",SocketsBuffer);
	
	// Get data format		
	WriteString(MySocket,"FORM:DATA?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);
	fprintf(fh,"Data form: %s",SocketsBuffer);	
	
	// Get system time before trace start
	WriteString(MySocket,"SYST:TIME?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);
	fprintf(fh,"Start time: %s",SocketsBuffer);

	// Initiate single sweep
	WriteString(MySocket,"INIT:CONT 0;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);
	printf("Trigger mode set? (1-Yes 0-No): %s",SocketsBuffer);

	WriteString(MySocket,"INIT:IMM;*OPC?\n");
	ReadString(MySocket,SocketsBuffer);
        //printf("%s",SocketsBuffer);

	for(i=0;i<atoi(argv[1]);i++)
	{
	  fallback:	WriteString(MySocket,"INIT:IMM;*OPC?\n");
				ReadString(MySocket,SocketsBuffer);
				WriteString(MySocket,"CALC:MEAS:WAOR?\n");
				ReadString(MySocket,SocketsBuffer);
				//printf("%s",SocketsBuffer);
				adc = atoi(SocketsBuffer);
				if( adc == 0 )
				{
					printf("Trace %d\n",++tcount);	
					WriteString(MySocket,"TRACE:DATA?\n");
					ReadTrace(MySocket,outfile,adc);
				}
				else if( adc == 1)
				{
					printf("Trace %d: ADC Over Range detected\n",++tcount);	
					WriteString(MySocket,"TRACE:DATA?\n");
					ReadTrace(MySocket,outfile,adc);
				}
	}

	// Get system time after trace ends
	WriteString(MySocket,"SYST:TIME?\n");
	ReadString(MySocket,SocketsBuffer);
	fprintf(fh,"End time: %s\n",SocketsBuffer);
	
SocketMainClose:
	
	// Close main port
	if(close(MySocket)==-1) 
	{
		printf("Error: Unable to close socket (%i)...\n",errno);
		perror("sockets"); // Print error message based on errno
		retval=1;;
	}

	exit(retval);
}
