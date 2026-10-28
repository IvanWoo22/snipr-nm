// Streaming FASTQ metrics and exact counts. C++17 + zlib; no external JSON dependency.
#include <zlib.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
namespace fs = std::filesystem;
using U = uint64_t;
using Hist = std::array<U,94>;
// Unlike gzgets(), this reader rejects trailing garbage after gzip members.
struct Fastq {
    std::ifstream file;
    z_stream z{};
    std::array<unsigned char,65536> input{}, output{};
    size_t position=0, available=0;
    bool compressed=false, source_eof=false, member_end=false;
    explicit Fastq(const std::string& name):file(name,std::ios::binary) {
        if(!file) throw std::runtime_error("Cannot open FASTQ: "+name);
        unsigned char magic[2]{};
        file.read(reinterpret_cast<char*>(magic),2);
        compressed=file.gcount()==2 && magic[0]==0x1f && magic[1]==0x8b;
        file.clear();file.seekg(0);
        if(!compressed && name.size()>=3 && name.substr(name.size()-3)==".gz")
            throw std::runtime_error("FASTQ .gz input has no gzip header");
        if(compressed && inflateInit2(&z,31)!=Z_OK) throw std::runtime_error("Cannot initialize gzip reader");
    }
    ~Fastq(){if(compressed)inflateEnd(&z);}
    size_t chunk() {
        if(!compressed){file.read(reinterpret_cast<char*>(output.data()),output.size());if(file.bad())throw std::runtime_error("FASTQ read error");return file.gcount();}
        for(;;){
            if(z.avail_in==0 && !source_eof){
                file.read(reinterpret_cast<char*>(input.data()),input.size());
                if(file.bad())throw std::runtime_error("FASTQ gzip read error");
                z.next_in=input.data();z.avail_in=static_cast<unsigned>(file.gcount());source_eof=file.eof();
            }
            if(member_end){
                if(z.avail_in==0 && source_eof)return 0;
                auto next=z.next_in;auto n=z.avail_in;
                if(inflateReset2(&z,31)!=Z_OK)throw std::runtime_error("Cannot reset gzip reader");
                z.next_in=next;z.avail_in=n;member_end=false;
            }
            z.next_out=output.data();z.avail_out=output.size();
            int status=inflate(&z,Z_NO_FLUSH);size_t produced=output.size()-z.avail_out;
            if(status==Z_STREAM_END)member_end=true;
            else if(status!=Z_OK && status!=Z_BUF_ERROR)throw std::runtime_error(std::string("Corrupt FASTQ gzip: ")+(z.msg?z.msg:"inflate error"));
            if(produced)return produced;
            if(source_eof && z.avail_in==0 && !member_end)throw std::runtime_error("Truncated FASTQ gzip stream");
        }
    }
    bool line(std::string& out){
        out.clear();
        for(;;){
            if(position==available){available=chunk();position=0;if(!available){if(out.empty())return false;break;}}
            const auto* begin=output.data()+position;
            const auto* newline=static_cast<const unsigned char*>(std::memchr(begin,'\n',available-position));
            const size_t length=newline?static_cast<size_t>(newline-begin):available-position;
            out.append(reinterpret_cast<const char*>(begin),length);
            position+=length;
            if(position<available){++position;break;}
        }
        if(!out.empty() && out.back()=='\r')out.pop_back();
        return true;
    }
};
struct Metrics {
    U total=0, bases=0, gc=0, qsum=0, q30=0, sampled=0, adapter_hits=0, low=0;
    U sample_limit;
    std::string adapter, adapter_label;
    std::map<size_t,U> lengths;
    std::vector<Hist> positions;
    std::unordered_set<std::string> sequences;
    Metrics(U limit, std::string a):sample_limit(limit),adapter(a),adapter_label(a) { for(char& c:adapter)c=std::toupper(static_cast<unsigned char>(c)); }
    static int code(char c) { switch(c) { case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;default:return -1; } }
    void add(std::string& seq, const std::string& qual, bool counting) {
        // Full-file validation and global statistics share one vectorizable pass.
        // ASCII case folding is valid only for ACGTN; reject all other bytes below.
        // Keep reductions local: no per-base locale calls, searches or shared writes.
        U read_gc=0, read_qsum=0, read_q30=0;
        unsigned invalid=0, lowercase=0;
        for(size_t i=0;i<seq.size();++i) {
            const unsigned original=static_cast<unsigned char>(seq[i]);
            const unsigned base=original&0xdfU;
            const unsigned quality=static_cast<unsigned char>(qual[i]);
            invalid |= !(base=='A'||base=='C'||base=='G'||base=='T'||base=='N');
            invalid |= quality<33 || quality>126;
            lowercase |= original>='a' && original<='z';
            read_gc += base=='C'||base=='G';
            read_qsum += quality-33;
            read_q30 += quality>=63;
        }
        if(counting && lowercase)throw std::runtime_error("Counter requires uppercase FASTQ sequence");
        if(invalid)throw std::runtime_error("Invalid sequence/quality at FASTQ record "+std::to_string(total+1));
        ++total; bases+=seq.size();++lengths[seq.size()];
        gc+=read_gc;qsum+=read_qsum;q30+=read_q30;
        if(total>sample_limit) return;
        // Only sampled sequence-dependent metrics need a normalized string.
        if(lowercase)for(char& c:seq)c=static_cast<char>(static_cast<unsigned char>(c)&0xdfU);
        ++sampled; sequences.insert(seq);
        adapter_hits+=!adapter.empty() && seq.find(adapter)!=std::string::npos;
        if(positions.size()<seq.size()) positions.resize(seq.size());
        for(size_t i=0;i<qual.size();++i) ++positions[i][qual[i]-33];
        std::array<U,16> dinuc{};U n=0;
        for(size_t i=1;i<seq.size();++i) {int a=code(seq[i-1]),b=code(seq[i]);if(a>=0&&b>=0){++dinuc[a*4+b];++n;}}
        double entropy=0; if(n) for(auto c:dinuc) if(c){double p=double(c)/n;entropy-=p*std::log2(p)/4;}
        low+=entropy<0.5;
    }
    static void ratio(std::ostream& o,U n,U d,double scale=100) {if(d)o<<scale*double(n)/d;else o<<"null";}
    void json(std::ostream& o) const {
        U inrange=0,mirna=0;for(auto [l,n]:lengths){if(l>=10&&l<=50)inrange+=n;if(l>=18&&l<=25)mirna+=n;}
        o<<std::setprecision(17)<<"{\"schema_version\":2,\"total_reads\":"<<total<<",\"total_bases\":"<<bases
         <<",\"reads_in_range\":"<<inrange<<",\"range_nt\":[10,50],\"median_length\":";
        if(!total)o<<"null";else{U at=0;double med=0;for(auto [l,n]:lengths){if(at<=(total-1)/2&&(total-1)/2<at+n)med+=l;if(at<=total/2&&total/2<at+n)med+=l;at+=n;}o<<med/2;}
        o<<",\"mean_quality\":";ratio(o,qsum,bases,1);o<<",\"q30_pct\":";ratio(o,q30,bases);
        o<<",\"gc_pct\":";ratio(o,gc,bases);o<<",\"mirna_18_25nt_pct\":";ratio(o,mirna,total);
        o<<",\"sampled_reads\":"<<sampled<<",\"sampling\":\"first N reads; not a random sample\",\"sequence_redundancy_pct\":";
        ratio(o,sampled-sequences.size(),sampled);o<<",\"exact_adapter_motif_pct\":";ratio(o,adapter_hits,sampled);
        o<<",\"adapter_motif\":\""<<adapter_label<<"\",\"low_dinucleotide_entropy_pct\":";ratio(o,low,sampled);
        o<<",\"low_entropy_definition\":\"Shannon entropy of ACGT dinucleotides / 4 < 0.5; N-containing pairs excluded\",\"length_distribution\":{";
        bool first=true;for(auto [l,n]:lengths){if(!first)o<<',';first=false;o<<'"'<<l<<"\":"<<n;}o<<"},\"quality_by_position\":{";
        for(size_t i=0;i<positions.size();++i){if(i)o<<',';U n=0,sum=0;for(unsigned q=0;q<94;++q){n+=positions[i][q];sum+=q*positions[i][q];}
            o<<'"'<<i+1<<"\":{\"reads\":"<<n<<",\"mean\":";ratio(o,sum,n,1);o<<",\"histogram\":{";first=true;
            for(unsigned q=0;q<94;++q)if(positions[i][q]){if(!first)o<<',';first=false;o<<'"'<<q<<"\":"<<positions[i][q];}o<<"}}";
        }o<<"}}\n";
    }
};
struct Counter {
    std::unordered_map<std::string,U> values;
    size_t limit,peak=0,created=0;U total=0,unique=0,spill_bytes=0;
    fs::path dir;std::vector<fs::path> chunks;
    explicit Counter(size_t n, fs::path p):limit(n),dir(std::move(p)) {}
    void spill() {
        if(values.empty())return;
        fs::create_directories(dir);
        auto path=dir/("chunk-"+std::to_string(created++)+".tsv");
        if(fs::exists(path))throw std::runtime_error("Refusing to overwrite count spill file");
        std::vector<const std::pair<const std::string,U>*> sorted;sorted.reserve(values.size());
        for(auto& x:values)sorted.push_back(&x);
        std::sort(sorted.begin(),sorted.end(),[](auto a,auto b){return a->first<b->first;});
        std::ofstream out(path);out.exceptions(std::ios::badbit|std::ios::failbit);
        for(auto p:sorted)out<<p->first<<'\t'<<p->second<<'\n';out.close();spill_bytes+=fs::file_size(path);chunks.push_back(path);values.clear();
    }
    void add(const std::string& seq){++values[seq];++total;peak=std::max(peak,values.size());if(values.size()>=limit)spill();}
    struct Item {std::string seq;U n;size_t i;bool operator>(const Item& b)const{return seq>b.seq;}};
    static bool next(std::ifstream& f,Item& item,size_t i){std::string line;if(!std::getline(f,line))return false;auto tab=line.find('\t');if(tab==std::string::npos)throw std::runtime_error("Bad count spill");item={line.substr(0,tab),std::stoull(line.substr(tab+1)),i};return true;}
    static U merge(const std::vector<fs::path>& paths,std::ostream& out) {
        std::vector<std::ifstream> inputs;for(auto& p:paths){inputs.emplace_back(p);if(!inputs.back())throw std::runtime_error("Missing count spill");inputs.back().exceptions(std::ios::badbit);}
        std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;
        for(size_t i=0;i<inputs.size();++i){Item x;if(next(inputs[i],x,i))queue.push(x);}
        std::string key;U count=0,unique=0;
        while(!queue.empty()){auto x=queue.top();queue.pop();if(!key.empty()&&key!=x.seq){out<<key<<'\t'<<count<<'\n';++unique;count=0;}key=x.seq;count+=x.n;Item y;if(next(inputs[x.i],y,x.i))queue.push(y);}
        if(!key.empty()){out<<key<<'\t'<<count<<'\n';++unique;}return unique;
    }
    void finish(const fs::path& output) {
        if(chunks.empty()) {
            std::vector<const std::pair<const std::string,U>*> sorted;sorted.reserve(values.size());for(auto& x:values)sorted.push_back(&x);
            std::sort(sorted.begin(),sorted.end(),[](auto a,auto b){return a->first<b->first;});
            std::ofstream out(output);out.exceptions(std::ios::badbit|std::ios::failbit);for(auto x:sorted)out<<x->first<<'\t'<<x->second<<'\n';unique=values.size();out.close();return;
        }
        spill();
        while(chunks.size()>32){std::vector<fs::path> reduced;
            for(size_t i=0;i<chunks.size();i+=32){std::vector<fs::path> group(chunks.begin()+i,chunks.begin()+std::min(i+32,chunks.size()));
                auto p=dir/("chunk-"+std::to_string(created++)+".tsv");if(fs::exists(p))throw std::runtime_error("Existing count spill");
                std::ofstream out(p);out.exceptions(std::ios::badbit|std::ios::failbit);merge(group,out);out.close();spill_bytes+=fs::file_size(p);reduced.push_back(p);for(auto& f:group)fs::remove(f);
            }chunks=std::move(reduced);
        }
        std::ofstream out(output);out.exceptions(std::ios::badbit|std::ios::failbit);unique=merge(chunks,out);out.close();for(auto& p:chunks)fs::remove(p);fs::remove(dir);
    }
};
int main(int argc,char** argv) {
    try {
        std::map<std::string,std::string> a;
        for(int i=1;i<argc;i+=2){if(i+1>=argc)throw std::runtime_error("Expected --option value pairs");a[argv[i]]=argv[i+1];}
        for(auto& [k,v]:a)if(k!="--input"&&k!="--metrics"&&k!="--adapter"&&k!="--sample-reads"&&k!="--counts"&&k!="--max-unique"&&k!="--spill-dir"&&k!="--count-meta")throw std::runtime_error("Unknown option: "+k);
        if(!a.count("--input")||!a.count("--metrics"))throw std::runtime_error("Required: --input FASTQ --metrics JSON");
        auto positive=[](const std::string& s){if(s.empty()||s.find_first_not_of("0123456789")!=std::string::npos||std::stoull(s)==0)throw std::runtime_error("Expected positive integer");return std::stoull(s);};
        U limit=positive(a.count("--sample-reads")?a["--sample-reads"]:"100000");
        U count_limit=positive(a.count("--max-unique")?a["--max-unique"]:"1000000");
        std::string adapter=a["--adapter"];for(char& c:adapter){c=std::toupper(static_cast<unsigned char>(c));if(std::string("ACGTN").find(c)==std::string::npos)throw std::runtime_error("Invalid adapter");}
        Metrics metrics(limit,a["--adapter"]);Counter counter(count_limit,a.count("--spill-dir")?a["--spill-dir"]:".count-spill");Fastq in(a["--input"]);
        std::string name,seq,plus,qual;bool count=a.count("--counts");
        while(in.line(name)) {
            if(!in.line(seq)||!in.line(plus)||!in.line(qual))throw std::runtime_error("Truncated FASTQ record "+std::to_string(metrics.total+1));
            if(name.empty()||name[0]!='@'||name.find_first_not_of(" \t\r\n",1)==std::string::npos||plus.empty()||plus[0]!='+'||seq.empty()||seq.size()!=qual.size())throw std::runtime_error("Malformed FASTQ record "+std::to_string(metrics.total+1));
            metrics.add(seq,qual,count);if(count)counter.add(seq);
        }
        if(count){auto output=a["--counts"]+".partial";counter.finish(output);fs::rename(output,a["--counts"]);}
        auto output=a["--metrics"]+".partial";std::ofstream out(output);out.exceptions(std::ios::badbit|std::ios::failbit);metrics.json(out);out.close();fs::rename(output,a["--metrics"]);
        if(a.count("--count-meta")){std::ofstream meta(a["--count-meta"]);meta.exceptions(std::ios::badbit|std::ios::failbit);meta<<"{\"total_reads\":"<<counter.total<<",\"unique_sequences\":"<<counter.unique<<",\"max_in_memory_unique\":"<<counter.peak<<",\"spill_files_created\":"<<counter.created<<",\"cumulative_spill_bytes\":"<<counter.spill_bytes<<"}\n";meta.close();}
        return 0;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
