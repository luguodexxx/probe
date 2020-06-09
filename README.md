## ProbeDesign: A RNA probe design tools [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

<pre>
__________              ___.          ________                .__               
\______   \_______  ____\_ |__   ____ \______ \   ____   _____|__| ____   ____  
 |     ___/\_  __ \/  _ \| __ \_/ __ \ |    |  \_/ __ \ /  ___/  |/ ___\ /    \ 
 |    |     |  | \(  <_> ) \_\ \  ___/ |    `   \  ___/ \___ \|  / /_/  >   |  \
 |____|     |__|   \____/|___  /\___  >_______  /\___  >____  >__\___  /|___|  /
                             \/     \/        \/     \/     \/  /_____/      \/ 
</pre>


ProbeDesign, a simple command line to generate probe. It's  ~~still not~~  completed.


### INSTALL

#### package install

python >= 3.6
```shell
git clone https://github.com/zhou-ran/probe
python setup.py install
```

#### third-party software

bowtie2

```shell
conda install bowtie2 -y

```

mfold
```shell

wget http://unafold.rna.albany.edu/download/mfold-3.6.tar.gz 
tar -xvf mfold-3.6.tar.gz
cd mfold-3.6
./configure --prefix=/path/to/save
make && make install
chmod 755 mfold
export PATH="/path/to/save/bin:$PATH"

# then add mfold to your env
```

[VarGibbs-2.2](http://bioinf.fisica.ufmg.br/software/vargibbs-2.2/)

### UASGE
There are three function.

#### transcript

Generate the probe based on the transcript ID.

```shell
# python ./scripts/trans_modify.py /path/to/cdna.fa /path/to/cdna.mod.fa
# bowtie2-bulk /path/to/cdna.mod.fa cDNAFastaFile

probedesign transcripts -faC cDNAFastaFile -config config.py
```

#### junction

Generate the probe based on the junction information.

```shell
probedesign junction -faG GenomeFastaFile -config config.py
```

#### circ

Generate the probe based on the circ-junction information.

```shell
probedesign circ -faG GenomeFastaFile -config config.py
```


### Reference
We generate the probe followed the same filter conditions with blockParse.py in [OligoMiner](https://github.com/brianbeliveau/OligoMiner)

