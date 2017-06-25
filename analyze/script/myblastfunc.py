import pandas as pd
from Bio.Blast import NCBIXML


def xml2df(filepath):
    """convert BLAST xml output to pd.Dataframe."""
    
    dct_lst=[]
    with open(filepath) as f:
        blastRecords = NCBIXML.parse(f)
        for rec in blastRecords:#one record for one query
            queryName=rec.query
            category=queryName.split('|')[0]
            queryLength=rec.query_length
            for alignment in rec.alignments:#one alignment record for one sbjct
                sbjctName=alignment.title
                sbjctLength=alignment.length
                for hsp in alignment.hsps:#one hsp for one (partial) alignment
                    dct={}
                    dct["category"]=category
                    dct["query_name"]=queryName
                    dct["sbjct_name"]=sbjctName
                    dct["evalue"]=hsp.expect
                    dct["query_length"]=queryLength
                    dct["sbjct_length"]=sbjctLength
                    dct["query_start"]=hsp.query_start
                    dct["quer_end"]=hsp.query_end
                    dct["sbjct_start"]=hsp.sbjct_start
                    dct["sbjct_end"]=hsp.sbjct_end
                    dct_lst.append(dct)
                    
    result_df=pd.DataFrame(dct_lst)
    result_df=result_df[["category","query_name","sbjct_name","evalue"]]
    return result_df


def tab2df(filepath):
    """convert BLAST tab output to pd.Dataframe.
    
        the naming of each column should be adjusted for the use in totalize_result method: category, evalue, query_name, etc...
    """

    #columns_lst=["qId", "tId", "seqIdentity", "alnLen", "mismatchCnt", "gapOpenCnt",
    #             "qStart", "qEnd", "tStart", "tEnd", "eVal", "bitScore"]
    #result_df=pd.read_csv(filepath, delimiter='\t', names=columns_lst)
    #result_df["category"]=[qId.split('|')[0] for qId in result_df["qId"]]
    
    columns_lst=["query_name","sbject_name","seqIdentity","alnLen","misMatchCnt","gapOpenCnt",
                 "query_start", "query_end", "sbjct_start","sbjct_end","evalue","bitScore"]
    result_df=pd.read_csv(filepath, delimiter='\t', names=columns_lst)
    result_df["category"]=[qId.split('|')[0][1:] for qId in result_df["query_name"]]#mmseqs has > in its name
    
    return result_df


def totalize_result(result_df,byQuery=False):
    """totalize #query for each category"""
    
    def format_count(count_df):
        """format count_df"""
        count_df=count_df.fillna(0)
        category_lst=get_category()
        count_df[category_lst]=count_df[category_lst].astype(int)
        count_df=count_df[["thres"]+category_lst]
        return count_df
    
    def get_category():
        category_lst=[]
        for prefix in ["F","DSF","DTSF"]:
            for i in range(2,7):
                category=prefix+str(i)
                category_lst.append(category)
        return category_lst
    
    thres_lst=[1,0.1,0.01,1e-3,1e-4,1e-5]
    dct_lst=[]
    for thres in thres_lst:
        dct={}
        dct["thres"]=thres
        category_lst=get_category()
        for category in category_lst:
            filtered_df=result_df[(result_df["evalue"]<thres) & (result_df["category"]==category)]
            if not(byQuery):
                dct[category]=filtered_df.shape[0]
            else:
                dct[category]=len(set(filtered_df["query_name"]))
        dct_lst.append(dct)
    count_df=pd.DataFrame(dct_lst)
    return format_count(count_df)