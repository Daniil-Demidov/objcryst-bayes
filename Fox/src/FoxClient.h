
#ifdef __WX__CRYST__
   // For compilers that support precompilation, includes "wx/wx.h".
   #ifndef __DARWIN__ // work around MacOSX type_info bug (??)
      #include "wx/wxprec.h"
   #endif

   #ifdef __BORLANDC__
       #pragma hdrstop
   #endif

   // for all others, include the necessary headers (this file is usually all you
   // need because it includes almost all "standard" wxWindows headers)
   #ifndef WX_PRECOMP
       #include "wx/wx.h"
   #endif

   //#include "wx/tooltip.h"
   //#include "wx/notebook.h"
   #include "wx/wfstream.h"
   #include "wx/zstream.h"
   #include "wx/fileconf.h"
   #include "wx/socket.h"
   #include "wx/process.h"
   #include "wx/sckstrm.h"
   #include "wx/sstream.h"
   #include "wx/zstream.h"
   #include "wx/wfstream.h"
   #include "wx/thread.h"
   #include "wx/stream.h"
   #include "wx/dir.h"
   #include "wx/stdpaths.h"
   //#include "wx/dynarray.h"
#endif

#include "wx/datetime.h"
#include "IOSocket.h"

#ifndef __FOX_CLIENT__
#define __FOX_CLIENT__

#define __CLIENT_LOGS 1
const wxEventType wxEVT_PROCESS_MY = wxNewEventType();

class ProcessMyEvent : public wxEvent
{
public:
    ProcessMyEvent(void* pSender)
    {
        SetId(-1);
        SetEventType(wxEVT_PROCESS_MY);
        m_sender = pSender;
        m_exit = false;
        m_pid = -1;
        m_status = -1;
    }

    virtual wxEvent* Clone() const
    {
        return new ProcessMyEvent(*this);
    }

    void*           m_sender;
    int             m_pid;
    int             m_status;
    wxString        m_dir;
    bool            m_exit;
};

struct ClientJob
{
   int ID;
   int nb_done;
   wxTimeSpan average_calc_time;
};

class FoxProcess
{
public:
    FoxProcess(wxString relDir);
    ~FoxProcess();

    void setPid(int pid);
    int  getPid();
    void     setTmpDir(wxString dir);
    wxString getTmpDir();
    void setRunning(bool run);
    bool isRunning();
    void setJobID(int id);
    int  getJobID();
    void setStarted(wxDateTime t);
    wxDateTime getStartingTime();
    int getProgressInPercents(wxTimeSpan avCalcTime);

private:

    int      pid;
    wxString tmpDIR;
    bool     running;
    int      jobID;
    wxDateTime startingtime;
};
class GrdRslt
{
public:
    GrdRslt(int ID, wxString cost, wxString content);
    ~GrdRslt();

    int ID;
    wxString cost;
    wxString content;
    bool sent;
	bool pending;
};



class FoxClient: public wxEvtHandler
{

public:
     FoxClient(wxString working_dir);
     ~FoxClient();
     bool ConnectClient(int nbOfTrial, wxString hostname);
     void OnTimerEvent(wxTimerEvent& event);
     void DoManyThingsOnTimer();
     //void OnSendResults(wxTimerEvent& event);
     void OnSocketEvent(wxSocketEvent &event);
     void WriteProtocol();
     bool IsClientConnected();
     void Disconnect();
     void OnProcessEvent(ProcessMyEvent& pEvent);
     void onProcessTerminate(int pid, int status, wxString dir);
     wxString getWorkingDir();

     //set nb of all available CPUs or cores on this PC
     void setNbOfAvailCPUs(int nb);

     //get nb of all available CPUs or cores on this PC
     //Do not use it for getting nb of unused CPUs!!!
     int  getNbOfAvailCPUs();

     //get nb of waiting processes (unused processes)
     //Do not use it for getting nb of all processes!!!
     int getNbOfUnusedProcesses();

     //Thread-safe way to get info about all processes
     //Returns copy of the processes, that can be used without mutex
     void get_copy_of_processes(vector<FoxProcess> &FP, vector<ClientJob> &CJ);

     //kill all running processes
     void KillProcesses();


     //bool   m_Connecting;
     bool   m_exit;
   void WriteMessageLog(wxString msg);

protected:
   wxString getJob(wxString inmsg, long pos);
   void SendCurrentState();
   void SaveResult(wxString fileName, wxString Cost, int ID, bool error);
   
   bool AnalyzeMessage(wxString msg);

   void SaveDataAsFile(wxString out, wxString filename);
   bool LoadFile(wxString filename, wxString &in);
   wxString getMyHostname();
   wxString addToPath(wxString str1, wxString str2);

   //reconnect client
   //void Reconnect();

   //it runs new job, return 0 if successful
   //job - content of the file, id - job id,
   int runNewJob(wxString job, int id, int nbTrial, bool rand);

   //send not-accepted jobs back to the server
   //void rejectJobs(std::vector<int> ids);

   //get file name "out-Cost-%f.xml", dir - directory with file
   wxString getOutputFile(wxString dir);

   //get cost from file name "out-Cost-%f.xml"
   wxString getCost(wxString filename);

   //reset processes (delete old and set new)
   void resetProcesses(int nbProcesses);

   //get only ONE of the unused processes
   FoxProcess *getUnusedProcess();

   wxSocketClient       * mpClient;
   wxString               m_hostname;
   wxTimer              * m_sendingTimer;
   vector<FoxProcess>     m_processes;
   vector<GrdRslt>        m_results;
   int                    m_nbOfAvailCPUs;
   IOSocket               m_IOSocket;
   wxString               m_working_dir;
   wxCriticalSection      m_ProcessCriticalSection;
   vector<ClientJob>      m_ListOfProcessedJobs;
   DECLARE_EVENT_TABLE()
};
class MyProcess : public wxProcess
{
public:
    MyProcess(FoxClient *parent, const wxString& cmd, wxString dir);
    virtual void OnTerminate(int pid, int status);

protected:
   wxString     m_cmd;
   FoxClient  * m_parent;
   wxString     m_dir;
};
#endif
