#! /usr/bin/env python3
import time
import re

def parse_speclist(iterable,size_limit=512000,recieve_timeout=120,encoding=None):
  data = []
  species_id = 0
  species_other_names_id = 0
  entry = dict()
  headerbarsre = re.compile(r'^=+$')
  headerspecre = re.compile(r'^\s*\(([0-9])\)\s+Real\s+organism\s+codes')
  headervirtre = re.compile(r'^\s*\(([0-9])\)\s+"Virtual"\s+codes')
  headersunderlinere = re.compile(r'^(\s*_+\s*)+$')
  nentryre = re.compile(r'^([A-Z0-9]{1,5})\s+([A-Z])\s+([0-9]+):\s+N=(.*)$')
  vnentryre = re.compile(r'^(9[A-Z0-9]{1,4})\s+([A-Z])\s+([0-9]+):\s+N=(.*)$')
  centryre = re.compile(r'^\s+C=(.*)$')
  sentryre = re.compile(r'^\s+S=(.*)$')
  i=0
  headernum = 0
  headeropen = False
  emptyrowcounter = 0
  headersread = False
  spec = None
  specset = set(('s','v'))
  parsedspec = set()
  rowcounter = 0
  size = 0
  start = time.time()
  chunkend = False
  remain = ''
  while True:
    try:
      
      chunk = next(iterable)
      if size_limit is not None:
        size += sys.getsizeof(chunk)
        if size > size_limit:
              raise StreamSizeLimitError('Response too large.')
      if recieve_timeout is not None and time.time() - start > recieve_timeout:
            raise StreamTimeoutError('Stream download time limit reached.')
      if encoding is not None:
        chunk = chunk.decode(encoding)
      chunk = remain+chunk
      
    except StopIteration:
      if chunkend:
        break
      else:
        lines = [remain]
        chunkend = True
        pass
    except:
      raise
    lines = chunk.split('\n')
    remain = lines.pop()
    
    for line in lines:
      
      if line != '':
        
        #spec code: 's' = read species, 'v' = read virtual taxons, 'o' read other (content ignored)
        #parsedspec code: 's' = parsed species, 'v' = parsed virtual taxons
        if headeropen:

          m1 = headerspecre.match(line)
          if m1:
              spec = 's'

              headernum = int(m1.group(1))

          else:
            m2 = headervirtre.match(line)
            if m2:
                spec = 'v'

                headernum = int(m2.group(1))

            elif headerbarsre.match(line) and spec is not None:

              headeropen = False
            elif spec is None:
                spec = 'o'
        else:
          if headerbarsre.match(line):
              headeropen = True


              if spec is not None:
                parsedspec.add(spec)

                spec = None
                headernum = 0
                
          elif spec == 's' and 's' not in parsedspec:
            if headernum == 1 and not headersread:

              if headersunderlinere.match(line):
                headersread = True

            else:
              nm = nentryre.match(line)
              if nm:

                if entry != dict():

                  data.append(entry)
                  entry = dict()

                  
                species_id += 1
                entry["id"] = species_id
                entry["code"] = nm.group(1)
                entry["kingdom"] = nm.group(2)
                entry["taxon_node"] = int(nm.group(3))
                entry["scientific_name"] = nm.group(4).strip()
                entry["common_name"] = []
                entry["synonym"] = []
              else:
                cm = centryre.match(line)
                if cm:
                  entry["common_name"].append(cm.group(1).strip())

                else:
                  sm = sentryre.match(line)
                  if sm:
                    entry["synonym"].append(sm.group(1).strip())

                  else:
                      #close section

                      data.append(entry)

                      parsedspec.add(spec)
                      spec = None
                      headernum = 0
          elif spec == 'v' and 'v' not in parsedspec:
            if headernum == 1 and not headersread:

              if headersunderlinere.match(line):

                headersread = True          
            else:
              nm = vnentryre.match(line)
              if nm:

                if entry != dict():

                  data.append(entry)
                  entry = dict()
                species_id += 1
                entry["id"] = species_id
                entry["code"] = nm.group(1)
                entry["kingdom"] = nm.group(2)
                entry["taxon_node"] = int(nm.group(3))
                entry["scientific_name"] = nm.group(4).strip()
                entry["common_name"] = []
                entry["synonym"] = []

              else:
                #close section

                data.append(entry)

                parsedspec.add(spec)
                spec = None
                headernum = 0
          elif specset.issubset(parsedspec):

            return data
       
  return data

def prepare_fixtures(data):
  data2 = {"uniprot_species":[],"uniprot_species_other_names":[]}
  
  entrys = {"model":"dynadb.DyndbUniprotSpecies", "pk":1,"fields":{"code":None,"kingdom":'O',"taxon_node":32630,"scientific_name":'synthetic'}}
  data2["uniprot_species"].append(entrys)
  entrys = {"model":"dynadb.DyndbUniprotSpecies", "pk":data[-1]["id"] + 1,"fields":{"code":None,"kingdom":'O',"taxon_node":32644,"scientific_name":'unidentified'}}
  data2["uniprot_species"].append(entrys)
  i = 1
  for entry in data:
    entrys = dict()
    entryo = dict()
    entrys["model"] = "dynadb.DyndbUniprotSpecies"
    entrys["pk"] = entry["id"] + 1
    entrys["fields"] = dict()
    entrys["fields"]["code"] = entry["code"]
    entrys["fields"]["kingdom"] = entry["kingdom"]
    entrys["fields"]["taxon_node"] =  entry["taxon_node"] 
    entrys["fields"]["scientific_name"] = entry["scientific_name"]
    if entrys != dict():
      data2["uniprot_species"].append(entrys)
    for name in entry["common_name"]:
      i += 1
      entryo["model"] = "dynadb.DyndbUniprotSpeciesAliases"
      entryo["pk"] = i
      entryo["fields"] = dict()
      entryo["fields"]["id_uniprot_species"] = entrys["pk"]
      entryo["fields"]["name"] = name
      entryo["fields"]["name_type"] = "C"
    if entryo != dict():
      data2["uniprot_species_other_names"].append(entryo)
      entryo = dict()
    for name in entry["synonym"]:
      entryo = dict()
      i += 1
      entryo["model"] = "dynadb.DyndbUniprotSpeciesAliases"
      entryo["pk"] = i
      entryo["fields"] = dict()
      entryo["fields"]["id_uniprot_species"] = entrys["pk"]
      entryo["fields"]["name"] = name
      entryo["fields"]["name_type"] = "S"
    if entryo != dict():
      data2["uniprot_species_other_names"].append(entryo)
  return data2
  

if __name__ == '__main__':
  import sys
  import os.path as path
  import json
  import numpy as np
  speclistfile = sys.argv[1]
  if len(sys.argv) > 2:
    outputpath = sys.argv[2]
  else:
    outputpath = './'
  if len(sys.argv) > 3:
    outputfilename = sys.argv[3]
  else:
    rootbasename,ext = path.splitext(path.basename(speclistfile))
  outputfilename = rootbasename + '.json'
  outputfilename1 = rootbasename + '_species.json'
  outputfilename2 = rootbasename + '_names.json'
  
  fout = open(path.join(outputpath,outputfilename),'w')
  fout1 = open(path.join(outputpath,outputfilename1),'w')
  fout2 = open(path.join(outputpath,outputfilename2),'w')
  fin = open(speclistfile,'r')
  data = parse_speclist(fin,size_limit=None,recieve_timeout=None)
  print('number of entries:'+str(len(data)))
  i=0
  for entry in data:
    i += len(entry['common_name'])
  print('number of common names:'+str(i))
  i=0
  for entry in data:
    i += len(entry['synonym'])
  print('number of synonym names:'+str(i))
  maxlen = np.max([np.max([len(entry["scientific_name"]),len(entry["common_name"]),len(entry["synonym"])]) for entry in data])
  print('Longest name has '+str(maxlen)+' characters.')
  fin.close()
  json.dump(data,fp=fout,indent=2)
  del data
  fout.close()
  fout = open(path.join(outputpath,outputfilename),'r')
  data = json.load(fp=fout)
  fout.close()
  data2 = prepare_fixtures(data)
  del data
  json.dump(data2["uniprot_species"],fp=fout1,indent=2)
  json.dump(data2["uniprot_species_other_names"],fp=fout2,indent=2)
  exit(0)
  