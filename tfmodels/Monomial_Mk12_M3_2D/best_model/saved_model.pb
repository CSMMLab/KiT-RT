خ<
��
W
AddN
inputs"T*N
sum"T"
Nint(0"!
Ttype:
2	��
D
AddV2
x"T
y"T
z"T"
Ttype:
2	��
B
AssignVariableOp
resource
value"dtype"
dtypetype�
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
~
BiasAddGrad
out_backprop"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
R
BroadcastGradientArgs
s0"T
s1"T
r0"T
r1"T"
Ttype0:
2	
N
Cast	
x"SrcT	
y"DstT"
SrcTtype"
DstTtype"
Truncatebool( 
T
CheckNumerics
tensor"T
output"T"
Ttype:
2"
messagestring�
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
8
Const
output"dtype"
valuetensor"
dtypetype
,
Exp
x"T
y"T"
Ttype:

2
^
Fill
dims"
index_type

value"T
output"T"	
Ttype"

index_typetype0:
2	
�
GatherV2
params"Tparams
indices"Tindices
axis"Taxis
output"Tparams"

batch_dimsint "
Tparamstype"
Tindicestype:
2	"
Taxistype:
2	
.
Identity

input"T
output"T"	
Ttype
,
Log
x"T
y"T"
Ttype:

2
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
>
Maximum
x"T
y"T
z"T"
Ttype:
2	
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(�
>
Minimum
x"T
y"T
z"T"
Ttype:
2	
?
Mul
x"T
y"T
z"T"
Ttype:
2	�
0
Neg
x"T
y"T"
Ttype:
2
	

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
�
Prod

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
@
ReadVariableOp
resource
value"dtype"
dtypetype�
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
?
Select
	condition

t"T
e"T
output"T"	
Ttype
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
H
ShardedFilename
basename	
shard

num_shards
filename
0
Sigmoid
x"T
y"T"
Ttype:

2
@
Softplus
features"T
activations"T"
Ttype:
2
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ��
@
StaticRegexFullMatch	
input

output
"
patternstring
�
StridedSlice

input"T
begin"Index
end"Index
strides"Index
output"T"	
Ttype"
Indextype:
2	"

begin_maskint "
end_maskint "
ellipsis_maskint "
new_axis_maskint "
shrink_axis_maskint 
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
<
Sub
x"T
y"T
z"T"
Ttype:
2	
�
Sum

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
P
	Transpose
x"T
perm"Tperm
y"T"	
Ttype"
Tpermtype0:
2	
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.7.02v2.7.0-rc1-69-gc256c071bb28��9
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
x
Adam/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/learning_rate
q
&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*
_output_shapes
: *
dtype0
h
VariableVarHandleOp*
_output_shapes
: *
dtype0*
shape:	*
shared_name
Variable
a
Variable/Read/ReadVariableOpReadVariableOpVariable*
_output_shapes
:	*
dtype0
p

Variable_1VarHandleOp*
_output_shapes
: *
dtype0*
shape
:		*
shared_name
Variable_1
i
Variable_1/Read/ReadVariableOpReadVariableOp
Variable_1*
_output_shapes

:		*
dtype0
�
layer_input/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:		�*#
shared_namelayer_input/kernel
z
&layer_input/kernel/Read/ReadVariableOpReadVariableOplayer_input/kernel*
_output_shapes
:		�*
dtype0
y
layer_input/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*!
shared_namelayer_input/bias
r
$layer_input/bias/Read/ReadVariableOpReadVariableOplayer_input/bias*
_output_shapes	
:�*
dtype0
�
block_0_layer_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameblock_0_layer_0/kernel
�
*block_0_layer_0/kernel/Read/ReadVariableOpReadVariableOpblock_0_layer_0/kernel* 
_output_shapes
:
��*
dtype0
�
block_0_layer_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameblock_0_layer_0/bias
z
(block_0_layer_0/bias/Read/ReadVariableOpReadVariableOpblock_0_layer_0/bias*
_output_shapes	
:�*
dtype0
�
block_1_layer_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameblock_1_layer_0/kernel
�
*block_1_layer_0/kernel/Read/ReadVariableOpReadVariableOpblock_1_layer_0/kernel* 
_output_shapes
:
��*
dtype0
�
block_1_layer_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameblock_1_layer_0/bias
z
(block_1_layer_0/bias/Read/ReadVariableOpReadVariableOpblock_1_layer_0/bias*
_output_shapes	
:�*
dtype0
�
block_2_layer_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameblock_2_layer_0/kernel
�
*block_2_layer_0/kernel/Read/ReadVariableOpReadVariableOpblock_2_layer_0/kernel* 
_output_shapes
:
��*
dtype0
�
block_2_layer_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameblock_2_layer_0/bias
z
(block_2_layer_0/bias/Read/ReadVariableOpReadVariableOpblock_2_layer_0/bias*
_output_shapes	
:�*
dtype0
�
block_3_layer_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameblock_3_layer_0/kernel
�
*block_3_layer_0/kernel/Read/ReadVariableOpReadVariableOpblock_3_layer_0/kernel* 
_output_shapes
:
��*
dtype0
�
block_3_layer_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameblock_3_layer_0/bias
z
(block_3_layer_0/bias/Read/ReadVariableOpReadVariableOpblock_3_layer_0/bias*
_output_shapes	
:�*
dtype0
�
block_4_layer_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameblock_4_layer_0/kernel
�
*block_4_layer_0/kernel/Read/ReadVariableOpReadVariableOpblock_4_layer_0/kernel* 
_output_shapes
:
��*
dtype0
�
block_4_layer_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameblock_4_layer_0/bias
z
(block_4_layer_0/bias/Read/ReadVariableOpReadVariableOpblock_4_layer_0/bias*
_output_shapes	
:�*
dtype0
�
block_5_layer_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameblock_5_layer_0/kernel
�
*block_5_layer_0/kernel/Read/ReadVariableOpReadVariableOpblock_5_layer_0/kernel* 
_output_shapes
:
��*
dtype0
�
block_5_layer_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameblock_5_layer_0/bias
z
(block_5_layer_0/bias/Read/ReadVariableOpReadVariableOpblock_5_layer_0/bias*
_output_shapes	
:�*
dtype0
�
dense_output/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*$
shared_namedense_output/kernel
|
'dense_output/kernel/Read/ReadVariableOpReadVariableOpdense_output/kernel*
_output_shapes
:	�*
dtype0
z
dense_output/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*"
shared_namedense_output/bias
s
%dense_output/bias/Read/ReadVariableOpReadVariableOpdense_output/bias*
_output_shapes
:*
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0
b
total_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_2
[
total_2/Read/ReadVariableOpReadVariableOptotal_2*
_output_shapes
: *
dtype0
b
count_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_2
[
count_2/Read/ReadVariableOpReadVariableOpcount_2*
_output_shapes
: *
dtype0
b
total_3VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_3
[
total_3/Read/ReadVariableOpReadVariableOptotal_3*
_output_shapes
: *
dtype0
b
count_3VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_3
[
count_3/Read/ReadVariableOpReadVariableOpcount_3*
_output_shapes
: *
dtype0
b
total_4VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_4
[
total_4/Read/ReadVariableOpReadVariableOptotal_4*
_output_shapes
: *
dtype0
b
count_4VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_4
[
count_4/Read/ReadVariableOpReadVariableOpcount_4*
_output_shapes
: *
dtype0
b
total_5VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_5
[
total_5/Read/ReadVariableOpReadVariableOptotal_5*
_output_shapes
: *
dtype0
b
count_5VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_5
[
count_5/Read/ReadVariableOpReadVariableOpcount_5*
_output_shapes
: *
dtype0
b
total_6VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_6
[
total_6/Read/ReadVariableOpReadVariableOptotal_6*
_output_shapes
: *
dtype0
b
count_6VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_6
[
count_6/Read/ReadVariableOpReadVariableOpcount_6*
_output_shapes
: *
dtype0
b
total_7VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_7
[
total_7/Read/ReadVariableOpReadVariableOptotal_7*
_output_shapes
: *
dtype0
b
count_7VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_7
[
count_7/Read/ReadVariableOpReadVariableOpcount_7*
_output_shapes
: *
dtype0
b
total_8VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_8
[
total_8/Read/ReadVariableOpReadVariableOptotal_8*
_output_shapes
: *
dtype0
b
count_8VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_8
[
count_8/Read/ReadVariableOpReadVariableOpcount_8*
_output_shapes
: *
dtype0
b
total_9VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_9
[
total_9/Read/ReadVariableOpReadVariableOptotal_9*
_output_shapes
: *
dtype0
b
count_9VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_9
[
count_9/Read/ReadVariableOpReadVariableOpcount_9*
_output_shapes
: *
dtype0
�
Adam/layer_input/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:		�**
shared_nameAdam/layer_input/kernel/m
�
-Adam/layer_input/kernel/m/Read/ReadVariableOpReadVariableOpAdam/layer_input/kernel/m*
_output_shapes
:		�*
dtype0
�
Adam/layer_input/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*(
shared_nameAdam/layer_input/bias/m
�
+Adam/layer_input/bias/m/Read/ReadVariableOpReadVariableOpAdam/layer_input/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/block_0_layer_0/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_0_layer_0/kernel/m
�
1Adam/block_0_layer_0/kernel/m/Read/ReadVariableOpReadVariableOpAdam/block_0_layer_0/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/block_0_layer_0/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_0_layer_0/bias/m
�
/Adam/block_0_layer_0/bias/m/Read/ReadVariableOpReadVariableOpAdam/block_0_layer_0/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/block_1_layer_0/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_1_layer_0/kernel/m
�
1Adam/block_1_layer_0/kernel/m/Read/ReadVariableOpReadVariableOpAdam/block_1_layer_0/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/block_1_layer_0/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_1_layer_0/bias/m
�
/Adam/block_1_layer_0/bias/m/Read/ReadVariableOpReadVariableOpAdam/block_1_layer_0/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/block_2_layer_0/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_2_layer_0/kernel/m
�
1Adam/block_2_layer_0/kernel/m/Read/ReadVariableOpReadVariableOpAdam/block_2_layer_0/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/block_2_layer_0/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_2_layer_0/bias/m
�
/Adam/block_2_layer_0/bias/m/Read/ReadVariableOpReadVariableOpAdam/block_2_layer_0/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/block_3_layer_0/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_3_layer_0/kernel/m
�
1Adam/block_3_layer_0/kernel/m/Read/ReadVariableOpReadVariableOpAdam/block_3_layer_0/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/block_3_layer_0/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_3_layer_0/bias/m
�
/Adam/block_3_layer_0/bias/m/Read/ReadVariableOpReadVariableOpAdam/block_3_layer_0/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/block_4_layer_0/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_4_layer_0/kernel/m
�
1Adam/block_4_layer_0/kernel/m/Read/ReadVariableOpReadVariableOpAdam/block_4_layer_0/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/block_4_layer_0/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_4_layer_0/bias/m
�
/Adam/block_4_layer_0/bias/m/Read/ReadVariableOpReadVariableOpAdam/block_4_layer_0/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/block_5_layer_0/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_5_layer_0/kernel/m
�
1Adam/block_5_layer_0/kernel/m/Read/ReadVariableOpReadVariableOpAdam/block_5_layer_0/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/block_5_layer_0/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_5_layer_0/bias/m
�
/Adam/block_5_layer_0/bias/m/Read/ReadVariableOpReadVariableOpAdam/block_5_layer_0/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/dense_output/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*+
shared_nameAdam/dense_output/kernel/m
�
.Adam/dense_output/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_output/kernel/m*
_output_shapes
:	�*
dtype0
�
Adam/dense_output/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*)
shared_nameAdam/dense_output/bias/m
�
,Adam/dense_output/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_output/bias/m*
_output_shapes
:*
dtype0
�
Adam/layer_input/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:		�**
shared_nameAdam/layer_input/kernel/v
�
-Adam/layer_input/kernel/v/Read/ReadVariableOpReadVariableOpAdam/layer_input/kernel/v*
_output_shapes
:		�*
dtype0
�
Adam/layer_input/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*(
shared_nameAdam/layer_input/bias/v
�
+Adam/layer_input/bias/v/Read/ReadVariableOpReadVariableOpAdam/layer_input/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/block_0_layer_0/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_0_layer_0/kernel/v
�
1Adam/block_0_layer_0/kernel/v/Read/ReadVariableOpReadVariableOpAdam/block_0_layer_0/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/block_0_layer_0/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_0_layer_0/bias/v
�
/Adam/block_0_layer_0/bias/v/Read/ReadVariableOpReadVariableOpAdam/block_0_layer_0/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/block_1_layer_0/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_1_layer_0/kernel/v
�
1Adam/block_1_layer_0/kernel/v/Read/ReadVariableOpReadVariableOpAdam/block_1_layer_0/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/block_1_layer_0/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_1_layer_0/bias/v
�
/Adam/block_1_layer_0/bias/v/Read/ReadVariableOpReadVariableOpAdam/block_1_layer_0/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/block_2_layer_0/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_2_layer_0/kernel/v
�
1Adam/block_2_layer_0/kernel/v/Read/ReadVariableOpReadVariableOpAdam/block_2_layer_0/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/block_2_layer_0/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_2_layer_0/bias/v
�
/Adam/block_2_layer_0/bias/v/Read/ReadVariableOpReadVariableOpAdam/block_2_layer_0/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/block_3_layer_0/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_3_layer_0/kernel/v
�
1Adam/block_3_layer_0/kernel/v/Read/ReadVariableOpReadVariableOpAdam/block_3_layer_0/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/block_3_layer_0/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_3_layer_0/bias/v
�
/Adam/block_3_layer_0/bias/v/Read/ReadVariableOpReadVariableOpAdam/block_3_layer_0/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/block_4_layer_0/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_4_layer_0/kernel/v
�
1Adam/block_4_layer_0/kernel/v/Read/ReadVariableOpReadVariableOpAdam/block_4_layer_0/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/block_4_layer_0/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_4_layer_0/bias/v
�
/Adam/block_4_layer_0/bias/v/Read/ReadVariableOpReadVariableOpAdam/block_4_layer_0/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/block_5_layer_0/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*.
shared_nameAdam/block_5_layer_0/kernel/v
�
1Adam/block_5_layer_0/kernel/v/Read/ReadVariableOpReadVariableOpAdam/block_5_layer_0/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/block_5_layer_0/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*,
shared_nameAdam/block_5_layer_0/bias/v
�
/Adam/block_5_layer_0/bias/v/Read/ReadVariableOpReadVariableOpAdam/block_5_layer_0/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/dense_output/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*+
shared_nameAdam/dense_output/kernel/v
�
.Adam/dense_output/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_output/kernel/v*
_output_shapes
:	�*
dtype0
�
Adam/dense_output/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*)
shared_nameAdam/dense_output/bias/v
�
,Adam/dense_output/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_output/bias/v*
_output_shapes
:*
dtype0
��
ConstConst*
_output_shapes
:	
�*
dtype0*��
value��B��	
�"��      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?      �?r��+s"�?�?���,�?`�}�?�YV��?{�]5v�? *��-�?��;�?q�?Q�?%�o�ކ�?%�o�ކ�?r�?Q�?��;�? *��-�?{�]5v�?�YV��?a�}�?�?���,�?t��+s"�?a��+s"���?���,��V�}���YV���{�]5v��
 *��-����;��r�?Q��%�o�ކ��%�o�ކ��r�?Q����;�� *��-��{�]5v���YV���k�}���?���,��}��+s"��2�e�;�?&���y�?>��*Ϳ?U�>n���?�a�A���?�-X���?����?"����+�?���8���?���8���?#����+�?����?�-X���?�a�A���?T�>n���?@��*Ϳ?"���y�?5�e�;�? �e�;�����y��3��*Ϳ�S�>n��ſ�a�A��ʿ�-X��ο���ѿ#����+ҿ���8��ҿ���8��ҿ#����+ҿ���ѿ�-X��ο�a�A��ʿZ�>n��ſJ��*Ϳ�%���y��?�e�;���r��M�?yv��5߽?0R�tc�?C!͌�?�U�<$g�?�+��Ǣ�?ܹ4�&�?z�f����?E6=|���?E6=|���?{�f����?ݹ4�&�?�+��Ǣ�?�U�<$g�?B!͌�?0R�tc�?sv��5߽?�r��M�?�r��M��lv��5߽��/R�tcȿA!͌п�U�<$gԿ�+��Ǣ׿ڹ4�&ڿ{�f���ۿE6=|��ܿE6=|��ܿ|�f���ۿ߹4�&ڿ�+��Ǣ׿�U�<$gԿF!͌п0R�tcȿwv��5߽��r��M����}h��?��l��?̊@d�?�P+����?���v���?�8�1�?����A�?qj?@`d�?��5����?��5����?rj?@`d�?����A�?�8�1�?���v���?�P+����?͊@d�?��l��?��}h��?��}h�����l�ÿƊ@dп�P+���տ���v��ڿ�8�1߿����A�rj?@`d���5������5����rj?@`d⿛���A῞8�1߿���v��ڿ�P+���տҊ@dп��l�ÿ��}h����@�,�?5�:u���?�o q���?��q(Ղ�?���iW�?r8��L��?�`)d���?��?Y�R�?�1�L��?�1�L��?��?Y�R�?�`)d���?r8��L��?���iW�?��q(Ղ�?�o q���?1�:u���?�@�,�?�@�,��+�:u��ǿo q��ӿ��q(Ղڿ���iW�u8��L�⿎`)d��俞�?Y�R濤1�L�翤1�L�翞�?Y�R濑`)d���s8��L�����iW���q(Ղڿ�o q��ӿ4�:u��ǿ�@�,������}�?�+��s�?5�(ói�?�X��1k�?T,��	��?�(�)���?�e�bE�?�+/����?둍�rj�?둍�rj�?�+/����?�e�bE�?�(�)���?T,��	��?�X��1k�?5�(ói�?�+��s�?����}�?����}����+��s˿,�(óiֿ�X��1k޿T,��	�⿈(�)����e�bE��+/����둍�rj�둍�rj��+/�����e�bE迅(�)���U,��	�⿧X��1k޿=�(óiֿ	�+��s˿����}����gy-U�?�'�0�?�fTA���??u[��?�S�ʞ�?�0a�?��?g(�0�m�?�����*�?���{��?���{��?����*�?h(�0�m�?�0a�?��?�S�ʞ�?>u[��?�fTA���?쓽'�0�?»gy-U�?��gy-U��哽'�0ο�fTA��ؿ>u[���S�ʞ��0a�?��e(�0�m�����*쿆��{������{�������*�j(�0�m��0a�?���S�ʞ�Cu[�࿖fTA��ؿ�'�0ο˻gy-U���ٿͫ��?v�����?���,�?�]�n<��?��X����?0��I ^�?hɌ��?"����?y*t<���?y*t<���?"����?iɌ��?0��I ^�?��X����?�]�n<��?���,�?s�����?�ٿͫ��?�ٿͫ���o����п���,ڿ�]�n<����X����4��I ^�fɌ��"�����y*t<���y*t<���"�����kɌ��1��I ^���X�����]�n<�����,ڿu����п�ٿͫ�����y�F;�?�:?�9��?;�Z@��?L����I�?q�sċ�?a0])V�?�j���?�&a��?�(�r��?�(�r��?�&a��?�j���?a0])V�?r�sċ�?K����I�? ;�Z@��?�:?�9��?��y�F;�?��y�F;���:?�9�п;�Z@�ڿJ����I�q�sċ�e0])V��j����&a���(�r���(�r���&a��k���b0])V�s�sċ�P����I�(;�Z@�ڿ�:?�9�п��y�F;��%�o�ކ�?r�?Q�?��;�? *��-�?{�]5v�?�YV��?`�}�?�?���,�?p��+s"�?g��+s"���?���,��^�}���YV���{�]5v�� *��-����;��r�?Q��%�o�ކ��%�o�ކ��r�?Q����;��	 *��-��{�]5v���YV���j�}���?���,��x��+s"��]��+s"�?�?���,�?U�}�?�YV��?{�]5v�? *��-�?��;�?r�?Q�?%�o�ކ�?���8���?#����+�?����?�-X���?�a�A���?V�>n���?>��*Ϳ?+���y�?0�e�;�?%�e�;��(���y��<��*Ϳ�V�>n��ſ�a�A��ʿ�-X��ο���ѿ#����+ҿ���8��ҿ���8��ҿ#����+ҿ���ѿ�-X��ο�a�A��ʿR�>n��ſI��*Ϳ�$���y��9�e�;���e�;�?���y�?1��*Ϳ?U�>n���?�a�A���?�-X���?����?#����+�?���8���?E6=|���?{�f����?ܹ4�&�?�+��Ǣ�?�U�<$g�?D!͌�?0R�tc�?�v��5߽?�r��M�?�r��M��|v��5߽�0R�tcȿD!͌п�U�<$gԿ�+��Ǣ׿ܹ4�&ڿ{�f���ۿE6=|��ܿE6=|��ܿ|�f���ۿ߹4�&ڿ�+��Ǣ׿�U�<$gԿ@!͌п0R�tcȿuv��5߽��r��M���r��M�?jv��5߽?�/R�tc�?C!͌�?�U�<$g�?+��Ǣ�?ڹ4�&�?{�f����?E6=|���?��5����?rj?@`d�?����A�?�8�1�?���v���?�P+����?̊@d�?l��?��}h��?��}h�����l�ÿˊ@dп�P+���տ���v��ڿ�8�1߿����A�rj?@`d���5������5����rj?@`d⿛���Aῡ8�1߿���v��ڿ�P+���տъ@dп��l�ÿ��}h�����}h��?��l��?Ɗ@d�?�P+����?���v���?�8�1�?����A�?rj?@`d�?��5����?�1�L��?��?Y�R�?�`)d���?s8��L��?���iW�?��q(Ղ�?�o q���?;�:u���?�@�,�?�@�,��8�:u��ǿ�o q��ӿ��q(Ղڿ���iW�s8��L�⿏`)d��俞�?Y�R濤1�L�翤1�L�翞�?Y�R濑`)d���t8��L�����iW࿽�q(Ղڿ�o q��ӿ2�:u��ǿ�@�,���@�,�?)�:u���?~o q���?��q(Ղ�?���iW�?q8��L��?�`)d���?��?Y�R�?�1�L��?둍�rj�?�+/����?�e�bE�?�(�)���?T,��	��?�X��1k�?5�(ói�?�+��s�?����}�?����}���+��s˿3�(óiֿ�X��1k޿T,��	�⿆(�)����e�bE��+/����둍�rj�둍�rj��+/�����e�bE过(�)���U,��	�⿛X��1k޿<�(óiֿ�+��s˿����}������}�?��+��s�?+�(ói�?�X��1k�?S,��	��?�(�)���?�e�bE�?�+/����?둍�rj�?���{��?����*�?g(�0�m�?�0a�?��?�S�ʞ�?@u[��?�fTA���?���'�0�?��gy-U�?��gy-U�����'�0ο�fTA��ؿ@u[���S�ʞ��0a�?��g(�0�m�����*쿆��{������{�������*�j(�0�m��0a�?���S�ʞ�=u[�࿕fTA��ؿ'�0οƻgy-U����gy-U�?㓽'�0�?�fTA���??u[��?�S�ʞ�?�0a�?��?e(�0�m�?����*�?���{��?y*t<���?"����?hɌ��?1��I ^�?��X����?�]�n<��?���,�?z�����?�ٿͫ��?�ٿͫ���x����п���,ڿ�]�n<����X����2��I ^�hɌ��"�����y*t<���y*t<���"�����kɌ��3��I ^���X�����]�n<�����,ڿt����п�ٿͫ����ٿͫ��?n�����?���,�?�]�n<��?��X����?/��I ^�?fɌ��?"����?y*t<���?�(�r��?�&a��?�j���?b0])V�?r�sċ�?M����I�?;�Z@��?�:?�9��?��y�F;�?��y�F;���:?�9�п;�Z@�ڿM����I�q�sċ�c0])V��j����&a���(�r���(�r���&a��k���d0])V�s�sċ�I����I�';�Z@�ڿ�:?�9�п��y�F;����y�F;�?�:?�9��?;�Z@��?L����I�?p�sċ�?`0])V�?�j���?�&a��?�(�r��?E�U��� ?�P?��oR?��* �h?��k �v?�>�`�3�?.�����?��˻:B�?uI砤�?��`.�?��`.�?wI砤�?��˻:B�?.�����?�>�`�3�?��k �v?��* �h?�P?��oR?H�U��� ?-�U��� ?�P?��oR?��* �h?��k �v?�>�`�3�?5�����?��˻:B�?wI砤�?��`.�?��`.�?wI砤�?��˻:B�?0�����?�>�`�3�?��k �v?��* �h?�P?��oR?U�U��� ?L/Ӆ�E? E��w?�,Ф��?�hOl�?�'f&Y�?��$C���?Ϗ`�+�?���
��?uɿV�?uɿV�?���
��?Ϗ`�+�?��$C���?�'f&Y�?�hOl�?�,Ф��? E��w?Q/Ӆ�E?//Ӆ�E? E��w?��,Ф��?�hOl�?�'f&Y�?�$C���?ˏ`�+�?���
��?uɿV�?uɿV�?���
��?я`�+�?��$C���?�'f&Y�?�hOl�?&�,Ф��? E��w?b/Ӆ�E?��a�KY?�}�d��?pM-d��?�0Rq�?=�O{s�?�h7�:u�?���o�^�?_IF�?�����?�����?`IF�?���o�^�?�h7�:u�?@�O{s�?}�0Rq�?rM-d��?�}�d��?��a�KY?��a�KY?�}�d��?bM-d��?{�0Rq�?=�O{s�?i7�:u�?���o�^�?`IF�?�����?�����?bIF�?���o�^�?�h7�:u�?B�O{s�?��0Rq�?~M-d��?�}�d��?-��a�KY?�I�:f?ʴ�`I�?N��M0�?r��l�ѽ?y7�N��?7����h�?gB1ќ�?,��EN$�?F�H��|�?F�H��|�?.��EN$�?gB1ќ�?7����h�?{7�N��?o��l�ѽ?P��M0�?´�`I�?�I�:f?�I�:f?���`I�?B��M0�?m��l�ѽ?y7�N��?A����h�?dB1ќ�?.��EN$�?F�H��|�?F�H��|�?.��EN$�?kB1ќ�?9����h�?}7�N��?���l�ѽ?Z��M0�?Ǵ�`I�?�I�:f?
�/[:p?�l#�u�?�eI	ٷ?��-K���?��3����?O��p
f�?��� k�?�gCe�$�?����<��?����<��?hCe�$�?��� k�?O��p
f�?��3����?��-K���?�eI	ٷ?�l#�u�?�/[:p?��/[:p?�l#�u�?�eI	ٷ?��-K���?��3����?V��p
f�?��� k�?hCe�$�?����<��?����<��?hCe�$�?��� k�?R��p
f�?��3����?��-K���?�eI	ٷ?�l#�u�?�/[:p?U�V"L]u?=�jQ�?:N��he�?��%�d��?��L��?��zE�|�?�J7j�?�K�H��?�F³\��?�F³\��?�K�H��?�J7j�?��zE�|�?��L��?��%�d��?:N��he�?4�jQ�?W�V"L]u?7�V"L]u?(�jQ�?!N��he�?��%�d��?��L��?��zE�|�?�J7j�?�K�H��?�F³\��?�F³\��?�K�H��?�J7j�?��zE�|�?��L��?��%�d��?PN��he�?9�jQ�?i�V"L]u?�.#��y?xc�^�{�?�meFT��?����S|�?��D�"��?�Pu���?WI����?~^X�i��?�G�yu_�?�G�yu_�?^X�i��?YI����?�Pu���?��D�"��?����S|�?�meFT��?mc�^�{�?�.#��y?�.#��y?`c�^�{�?rmeFT��?����S|�?��D�"��?�Pu���?TI����?^X�i��?�G�yu_�?�G�yu_�?�^X�i��?\I����?�Pu���?��D�"��?���S|�?�meFT��?uc�^�{�?(�.#��y?����#}?ֻ+��?NҖ�=i�?o���.��?�
�>��?��.Q'�?`�d5��?1-�JM��?!x�.���?!x�.���?2-�JM��?b�d5��?��.Q'�?��
�>��?m���.��?OҖ�=i�?ֻ+��?����#}?����#}?ֻ+��??Җ�=i�?k���.��?�
�>��?��.Q'�?\�d5��?2-�JM��?!x�.���?!x�.���?4-�JM��?e�d5��?��.Q'�?��
�>��?x���.��?^Җ�=i�?ֻ+��?���#}?�I�|��~?�h��?�d����?~���n��?"a�. ��?����hQ�?=AB��?��AO��?���qX��?���qX��?��AO��??AB��?����hQ�?%a�. ��?{���n��?�d����?��h��?�I�|��~?�I�|��~?�h��?�d����?y���n��?"a�. ��?����hQ�?9AB��?��AO��?���qX��?���qX��?��AO��?BAB��?����hQ�?(a�. ��?����n��?�d����?�h��?�I�|��~?XS�u�W?�>�`�3q?�����Zz?&�*�?�>�`�3�?&�*�?�����Zz?�>�`�3q?VS�u�W?MS�u�W��>�`�3q������Zz�&�*���>�`�3��&�*�������Zz�}>�`�3q�ZS�u�W�GS�u�W?y>�`�3q?�����Zz?&�*�?�>�`�3�?&�*�?�����Zz?~>�`�3q?^S�u�W?CS�u�W�x>�`�3q������Zz�&�*���>�`�3��&�*�������Zz�>�`�3q�dS�u�W��Y�7�~?�'f&Y�?@��C��?|���Ȥ?�'f&Y�?}���Ȥ?@��C��?�'f&Y�?�Y�7�~?�Y�7�~��'f&Y��>��C��}���Ȥ��'f&Y��|���Ȥ�A��C�񠿬'f&Y���Y�7�~��Y�7�~?�'f&Y�?;��C��?|���Ȥ?�'f&Y�?|���Ȥ?C��C��?�'f&Y�?�Y�7�~?�Y�7�~��'f&Y��:��C��|���Ȥ��'f&Y������Ȥ�D��C�񠿯'f&Y���Y�7�~�-/ M�??�O{s�?�f�3�?w-��r�??�O{s�?x-��r�?�f�3�?D�O{s�?+/ M�?$/ M��A�O{s���f�3x-��r��?�O{s��w-��r���f�3:�O{s��./ M��/ M�?4�O{s�?�f�3�?v-��r�?@�O{s�?v-��r�?�f�3�?;�O{s�?1/ M�?/ M��3�O{s���f�3w-��r��?�O{s��z-��r���f�3=�O{s��6/ M��b����z�?y7�N䨶?����[�?�nT$K�?z7�N��?�nT$K�?����[�?~7�N䨶?_����z�?S����z��}7�N䨶�����[���nT$Kſz7�N�ƿ�nT$Kſ����[��v7�N䨶�e����z��L����z�?p7�N䨶?����[�?�nT$K�?{7�N��?�nT$K�?����[�?w7�N䨶?k����z�?F����z��o7�N䨶�����[���nT$Kſy7�N�ƿ�nT$Kſ����[��x7�N䨶�r����z��k�E�/�?��3����?J�r%��?��}�	^�?��3����?��}�	^�?J�r%��?��3����?i�E�/�?_�E�/����3�����H�r%�ɿ��}�	^Ͽ��3���п��}�	^ϿK�r%�ɿ��3�����m�E�/��[�E�/�?��3����?C�r%��?��}�	^�?��3����?��}�	^�?P�r%��?��3����?p�E�/�?V�E�/����3�����B�r%�ɿ��}�	^Ͽ��3���п��}�	^ϿR�r%�ɿ��3�����v�E�/���V�^K��?��L��?�8���?�.Dۥ�?��L��?�.Dۥ�?�8���?��L��?�V�^K��?�V�^K�����L�ſ�8��п�.DۥԿ��L�տ�.DۥԿ�8��п��L�ſ�V�^K����V�^K��?��L��?�8���?�.Dۥ�?��L��?�.Dۥ�?�8���?��L��?�V�^K��?�V�^K�����L�ſ�8��п�.DۥԿ��L�տ�.DۥԿ�8��п��L�ſ�V�^K���Zf-�hu�?��D�"��?[u�:�[�?��{f���?��D�"��?��{f���?[u�:�[�?��D�"��?Xf-�hu�?Pf-�hu����D�"�ʿZu�:�[Կ��{f��ؿ��D�"�ڿ��{f��ؿ[u�:�[Կ��D�"�ʿ[f-�hu��Lf-�hu�?��D�"��?Uu�:�[�?��{f���?��D�"��?��{f���?`u�:�[�?��D�"��?^f-�hu�?If-�hu����D�"�ʿUu�:�[Կ��{f��ؿ��D�"�ڿ��{f��ؿau�:�[Կ��D�"�ʿcf-�hu��/�e�Ѵ?��
�>��?�=�G��?d�H�)�?��
�>��?d�H�)�?�=�G��?��
�>��?-�e�Ѵ?&�e�Ѵ���
�>�Ϳ�=�G�ֿd�H�)ܿ��
�>�ݿd�H�)ܿ�=�G�ֿ�
�>�Ϳ1�e�Ѵ�!�e�Ѵ?�
�>��?�=�G��?d�H�)�?��
�>��?d�H�)�?�=�G��?�
�>��?5�e�Ѵ?�e�Ѵ��
�>�Ϳ�=�G�ֿd�H�)ܿ��
�>�ݿd�H�)ܿ�=�G�ֿ�
�>�Ϳ:�e�Ѵ�v���5�?$a�. ��?U��UV�?2�	b���?#a�. ��?3�	b���?U��UV�?*a�. ��?t���5�?k���5��'a�. �ϿT��UVؿ3�	b��ݿ#a�. �߿2�	b��ݿV��UVؿa�. �Ͽx���5��f���5�?a�. ��?O��UV�?1�	b���?%a�. ��?1�	b���?[��UV�? a�. ��?|���5�?b���5��a�. �ϿN��UVؿ2�	b��ݿ#a�. �߿6�	b��ݿ\��UVؿ"a�. �Ͽ����5����`.�?wI砤�?��˻:B�?0�����?�>�`�3�?��k �v?��* �h?�P?��oR?C�U��� ?6�U��� ?�P?��oR?��* �h?��k �v?�>�`�3�?2�����?��˻:B�?wI砤�?��`.�?��`.�?wI砤�?��˻:B�?3�����?�>�`�3�?��k �v?��* �h?�P?��oR?N�U��� ?'�U��� ?�P?��oR?��* �h?��k �v?~>�`�3�?-�����?��˻:B�?wI砤�?��`.�?uɿV�?���
��?Ϗ`�+�?��$C���?�'f&Y�?�hOl�?�,Ф��?" E��w?I/Ӆ�E?7/Ӆ�E? E��w?�,Ф��?�hOl�?�'f&Y�?�$C���?Ϗ`�+�?���
��?uɿV�?uɿV�?���
��?я`�+�?�$C���?�'f&Y�?�hOl�?$�,Ф��? E��w?X/Ӆ�E?'/Ӆ�E?��D��w?��,Ф��?�hOl�?�'f&Y�?��$C���?ˏ`�+�?���
��?uɿV�?�����?`IF�?���o�^�?�h7�:u�?@�O{s�?��0Rq�?pM-d��?�}�d��?��a�KY?���a�KY?�}�d��?mM-d��?��0Rq�?=�O{s�?i7�:u�?���o�^�?`IF�?�����?�����?bIF�?���o�^�?i7�:u�?B�O{s�?y�0Rq�?|M-d��?�}�d��? ��a�KY?��a�KY?�}�d��?aM-d��?�0Rq�?;�O{s�?�h7�:u�?���o�^�?`IF�?�����?F�H��|�?.��EN$�?gB1ќ�?9����h�?{7�N��?u��l�ѽ?N��M0�?ִ�`I�?�I�:f?�I�:f?Ѵ�`I�?L��M0�?u��l�ѽ?y7�N��?=����h�?gB1ќ�?.��EN$�?F�H��|�?F�H��|�?.��EN$�?kB1ќ�??����h�?}7�N��?j��l�ѽ?X��M0�?Ŵ�`I�?�I�:f?�I�:f?���`I�?B��M0�?r��l�ѽ?v7�N��?5����h�?dB1ќ�?.��EN$�?F�H��|�?����<��?hCe�$�?��� k�?R��p
f�?��3����?��-K���?�eI	ٷ?�l#�u�?�/[:p?��/[:p?�l#�u�?�eI	ٷ?��-K���?��3����?R��p
f�?��� k�?hCe�$�?����<��?����<��?hCe�$�?��� k�?T��p
f�?��3����?��-K���?�eI	ٷ?�l#�u�?�/[:p?�/[:p?�l#�u�?�eI	ٷ?��-K���?��3����?M��p
f�?��� k�?hCe�$�?����<��?�F³\��?�K�H��?�J7j�?��zE�|�?��L��?��%�d��?:N��he�?I�jQ�?P�V"L]u?@�V"L]u?B�jQ�?4N��he�?��%�d��?��L��?��zE�|�?�J7j�?�K�H��?�F³\��?�F³\��?�K�H��?�J7j�?��zE�|�?��L��?~�%�d��?NN��he�?8�jQ�?`�V"L]u?-�V"L]u?&�jQ�?N��he�?��%�d��?��L��?��zE�|�?�J7j�?�K�H��?�F³\��?�G�yu_�?^X�i��?WI����?�Pu���?��D�"��? ���S|�?�meFT��?�c�^�{�?	�.#��y?�.#��y?~c�^�{�?~meFT��? ���S|�?��D�"��?�Pu���?WI����?^X�i��?�G�yu_�?�G�yu_�?�^X�i��?\I����?�Pu���?��D�"��?����S|�?�meFT��?qc�^�{�?�.#��y?ށ.#��y?\c�^�{�?rmeFT��?����S|�?��D�"��?�Pu���?TI����?^X�i��?�G�yu_�?!x�.���?2-�JM��?`�d5��?��.Q'�?��
�>��?q���.��?NҖ�=i�?ֻ+��?����#}?����#}?ֻ+��?KҖ�=i�?q���.��?�
�>��?��.Q'�?`�d5��?2-�JM��?!x�.���?!x�.���?4-�JM��?e�d5��?��.Q'�?��
�>��?i���.��?\Җ�=i�?ֻ+��?���#}?����#}?ֻ+��?=Җ�=i�?o���.��?�
�>��?��.Q'�?\�d5��?2-�JM��?!x�.���?���qX��?��AO��?=AB��?����hQ�?%a�. ��?����n��?�d����?��h��?�I�|��~?�I�|��~?��h��?�d����?����n��?"a�. ��?����hQ�?=AB��?��AO��?���qX��?���qX��?��AO��?BAB��?����hQ�?(a�. ��?w���n��?�d����?�h��?�I�|��~?�I�|��~?�h��?�d����?~���n��?a�. ��?����hQ�?9AB��?��AO��?���qX��?�O�<�.�>4u��v�?nB�o�%?3�`3��:?�e� N9I?lobA~�S?�(j:1�Z?�}��`?{o���a?{o���a?�}��`?�(j:1�Z?lobA~�S?�e� N9I?/�`3��:?pB�o�%?*u��v�?�O�<�.�>tO�<�.��u��v��WB�o�%�+�`3��:��e� N9I�uobA~�S��(j:1�Z��}��`�{o���a�{o���a��}��`��(j:1�Z�nobA~�S��e� N9I�H�`3��:��B�o�%�1u��v���O�<�.��=
��c��>0(\MC�<?\Y�oh_?O�{�%�s?s'>��c�?�۸�Ɩ�?٣�U\�?�!�8�o�?�H�׵�?�H�׵�?�!�8�o�?٣�U\�?�۸�Ɩ�?u'>��c�?L�{�%�s?bY�oh_?(\MC�<?C
��c��>
��c��(\MC�<�;Y�oh_�I�{�%�s�s'>��c���۸�Ɩ��ң�U\���!�8�o���H�׵���H�׵���!�8�o��ܣ�U\���۸�Ɩ��w'>��c��]�{�%�s�Y�oh_�+(\MC�<�W
��c��g�`��?�i�Z?_�e U|?\f0+��?E[����?��0!ʩ?b��v�?>䝙l$�?wB�l[1�?wB�l[1�?@䝙l$�?d��v�?��0!ʩ?H[����?Yf0+��?_�e U|?�i�Z?�g�`��?Rg�`����i�Z��^�e U|�Vf0+���E[������0!ʩ�^��v��@䝙l$��wB�l[1��wB�l[1��B䝙l$��h��v����0!ʩ�J[�����ff0+���$_�e U|��i�Z��g�`����j(�H"?� f[	�m?M��H�?�"�4�Z�?{&^gn�?�MD��?��cd �?H�&zM�?FT���?FT���? H�&zM�?��cd �?�MD��?~&^gn�?�"�4�Z�?P��H�?� f[	�m?�j(�H"?hj(�H"�� f[	�m�;��H���"�4�Z��{&^gn���MD�����cd Ŀ H�&zMȿFT��ʿFT��ʿ H�&zMȿ��cd Ŀ�MD����&^gn���"�4�Z��`��H��� f[	�m��j(�H"�D��f�W0?� ���z?#L�{]�??6%�<2�?.?����?�sKce��?�<g@��?�br�ع�?�u6aE��?�u6aE��?�br�ع�?�<g@��?�sKce��?1?����?;6%�<2�?(L�{]�?� ���z?J��f�W0?#��f�W0�} ���z�L�{]��96%�<2��.?������sKce�ʿ�<g@�ѿ�br�عտ�u6aE�׿�u6aE�׿�br�عտ�<g@�ѿ�sKce�ʿ1?�����I6%�<2��CL�{]��� ���z�]��f�W0�h��� �8?{|	�4�?*�m��?�J2�|�?�����?�E�)�?��N���? �=��h�?x,#%5 �?x,#%5 �?"�=��h�?��N���?�E�)�?�����?�J2�|�?*�m��?p|	�4�?l��� �8?4��� �8�`|	�4���)�m����J2�|������ɿ��E�)Կ��N��ۿ"�=��h�x,#%5 �x,#%5 �$�=��h���N��ۿ�E�)Կ����ɿ�J2�|��.*�m���v|	�4������ �8���k��j@?u�uB1ߊ?�"؂q?�?B?��G�?����?Ȫ{fj��?�x0�A�?�N��F��?o��*��?o��*��?�N��F��?�x0�A�?Ȫ{fj��?�����???��G�?�"؂q?�?e�uB1ߊ?��k��j@?\�k��j@�R�uB1ߊ��"؂q?��??��G¿���ѿ֪{fj�ڿ�x0�A⿾N��F��o��*��o��*���N��F�志x0�A�̪{fj�ڿ����ѿP?��G¿#؂q?��p�uB1ߊ���k��j@��g����C?PjV�?2�S����?����P��?g���,��?̳�jK��?3e��U��?a0rr#�?Ha��m��?Ha��m��?a0rr#�?6e��U��?̳�jK��?i���,��?����P��?4�S����?�PjV�?�g����C?�g����C��PjV�� �S�����~���P�ſg���,�Կ۳�jK�߿/e��U��a0rr#�Ha��m��Ha��m��a0rr#�:e��U��г�jK�߿l���,�Կ����P�ſF�S�����PjV��h����C��gY��uE?hD�	��?5����?`�[����?1Ib�?&�-�Zf�?��ʐ�?�q�[��?7�Ȯ�K�?7�Ȯ�K�?�q�[��?��ʐ�?&�-�Zf�?4Ib�?\�[����?7����?uhD�	��?�gY��uE?sgY��uE�ihD�	��������X�[���ǿ1Ibֿ.�-�Zf���ʐ��q�[��7�Ȯ�K�7�Ȯ�K��q�[����ʐ�(�-�Zf�7Ibֿo�[���ǿH�����|hD�	����gY��uE�mJ�H�F�>�Rw"?��(l�7?��G:C?�e� N9I?�
�_uK?�j�'�H?���k:A?H:���(?>:���(����k:A��j�'�H��
�_uK��e� N9I���G:C���(l�7�	�Rw"�pJ�H�F�TJ�H�F� �Rw"�r�(l�7���G:C��e� N9I��
�_uK��j�'�H����k:A�Q:���(�4:���(?���k:A?�j�'�H?�
�_uK?�e� N9I? �G:C?��(l�7?�Rw"?~J�H�F�>6���0)?��o�h�Z?A�Ps��p?S3Kf�|?s'>��c�?I��w��?�v��'�? �ȴy?�53���a?�53���a� �ȴy��v��'��I��w���t'>��c��Q3Kf�|�B�Ps��p���o�h�Z�<���0)����0)���o�h�Z�6�Ps��p�O3Kf�|�t'>��c��K��w����v��'�� �ȴy��53���a��53���a? �ȴy?�v��'�?I��w��?t'>��c�?^3Kf�|?K�Ps��p?��o�h�Z?P���0)?oP5�F?/�da�Ix?<��0a�?`�"��I�?F[����?0 g��?#b9�I�?�.����?A����;�?:����;���.�����#b9�I��0 g���G[�����^�"��I��?��0a��%�da�Ix�sP5�F�LP5�F��da�Ix�)��0a��\�"��I��G[�����0 g���#b9�I���.�����G����;��4����;�?�.����?�"b9�I�?0 g��?H[����?h�"��I�?Q��0a�?,�da�Ix?�P5�F?ar���Z?`�pQ��?O��j�u�?=����?|&^gn�?�������?�i�Ӹ�?YS��(�?0=��˨�?)=��˨��XS��(���i�Ӹ����������}&^gn��<�����Q��j�u��W�pQ�ꋿfr���Z�<r���Z�I�pQ�ꋿD��j�u��;�����|&^gn�����������i�Ӹ��RS��(��7=��˨��"=��˨�?HS��(�?�i�Ӹ�?�������?}&^gn�?H����?[��j�u�?]�pQ��?{r���Z?h�dD�Yg?�i�R)��?3���7�?������?/?����?n�9|��?���\���?��3+I�?B�V��?B�V�����3+I�����\����n�9|�¿0?������������5���7���i�R)���n�dD�Yg�H�dD�Yg��i�R)������7���������/?�����o�9|�¿���\������3+I��B�V���B�V��?��3+I�?���\���?p�9|��?.?����?������?G���7�?�i�R)��?�dD�Yg?�0`�ߢq?�u��٢?@,p�+��?S�UW���?�����?���g��?�'�KNH�?W��c��?|�꟧2�?r�꟧2��U��c����'�KNHɿ���g�̿����ɿR�UW��ÿ@,p�+����u��٢��0`�ߢq�}0`�ߢq��u��٢�/,p�+���R�UW��ÿ����ɿ���g�̿�'�KNHɿQ��c�����꟧2��h�꟧2�?M��c��?�'�KNH�?���g��?�����?[�UW���?N,p�+��?�u��٢?�0`�ߢq?�aauw?J�]�_�?�@��o\�?I�Ɓ�?�����?���4��?Ɨ}B��?�&6�ld�?���v���?���v�����&6�ldǿƗ}B�п���4�ҿ����ѿH�Ɓʿ�@��o\��@�]�_���aauw��aauw�5�]�_���@��o\��I�Ɓʿ����ѿ���4�ҿƗ}B�п�&6�ldǿ��v�����v���?�&6�ld�?�ŗ}B��?���4��?�����?S�Ɓ�?�@��o\�?G�]�_�?bauw?�o,s�|?YO���?������?��e��C�?h���,��?��t�LS�?��_��"�?�V��?�WMK��?�WMK����V�̿��_��"Կ��t�LSֿh���,�Կ��e��CϿ�����¿ YO�����o,s�|�oo,s�|��XO���������¿��e��CϿi���,�Կ��t�LSֿ��_��"Կ�V�̿�WMK����WMK��?�V��?��_��"�?��t�LS�?i���,��?��e��C�?������?YO���?�o,s�|?|�ck�~?�n�/�b�?�|�d&�?�]gK��?2Ib�?�/�]�?��SF��?�Pl,���?b�x+t�?Y�x+t絿�Pl,��ο��SF�տ�/�]ؿ3Ibֿ�]gK�ѿ�|�d&Ŀ�n�/�b����ck�~�O�ck�~��n�/�b���|�d&Ŀ�]gK�ѿ3Ibֿ�/�]ؿ��SF�տ�Pl,��οj�x+t絿P�x+t�?�Pl,���?���SF��?�/�]�?4Ib�?�]gK��?�|�d&�?�n�/�b�?��ck�~?J:���(?���k:A?�j�'�H?�
�_uK?�e� N9I?��G:C?��(l�7?�Rw"?kJ�H�F�>^J�H�F�>�Rw"?}�(l�7?��G:C?�e� N9I?�
�_uK?�j�'�H?���k:A?L:���(?8:���(����k:A��j�'�H��
�_uK��e� N9I���G:C���(l�7��Rw"�wJ�H�F�NJ�H�F���Rw"�p�(l�7���G:C��e� N9I��
�_uK��j�'�H����k:A�V:���(��53���a? �ȴy?�v��'�?I��w��?t'>��c�?T3Kf�|?A�Ps��p?��o�h�Z?3���0)?���0)?��o�h�Z??�Ps��p?T3Kf�|?s'>��c�?J��w��?�v��'�? �ȴy?�53���a?�53���a� �ȴy��v��'��I��w���u'>��c��N3Kf�|�J�Ps��p���o�h�Z�D���0)����0)���o�h�Z�4�Ps��p�S3Kf�|�q'>��c��K��w����v��'�� �ȴy��53���a�B����;�?�.����?#b9�I�?0 g��?G[����?c�"��I�?<��0a�?:�da�Ix?jP5�F?WP5�F?4�da�Ix?8��0a�?c�"��I�?F[����?0 g��?#b9�I�?�.����?D����;�?6����;���.����� #b9�I��0 g���H[�����Z�"��I��M��0a��(�da�Ix�yP5�F�EP5�F��da�Ix�'��0a��`�"��I��F[�����0 g���#b9�I���.�����K����;��2=��˨�?TS��(�?�i�Ӹ�?�������?}&^gn�??����?O��j�u�?m�pQ��?[r���Z?Hr���Z?h�pQ��?M��j�u�??����?|&^gn�?�������?�i�Ӹ�?PS��(�?4=��˨�?%=��˨��JS��(���i�Ӹ����������~&^gn��9�����Y��j�u��[�pQ�ꋿor���Z�3r���Z�F�pQ�ꋿD��j�u��=�����z&^gn�����������i�Ӹ��SS��(��<=��˨��B�V��?��3+I�?���\���?p�9|��?0?����?������?3���7�?�i�R)��?e�dD�Yg?Q�dD�Yg?�i�R)��?.���7�?������?/?����?n�9|��?���\���?��3+I�?B�V��?B�V�����3+I�����\����o�9|�¿0?������������C���7���i�R)���t�dD�Yg�@�dD�Yg��i�R)������7���������,?�����o�9|�¿���\������3+I��"B�V����꟧2�?S��c��?�'�KNH�?���g��?�����?U�UW���?@,p�+��?$�u��٢?�0`�ߢq?�0`�ߢq?�u��٢?;,p�+��?U�UW���?�����?���g��?�'�KNH�?P��c��?��꟧2�?m�꟧2��M��c����'�KNHɿ���g�̿����ɿP�UW��ÿM,p�+����u��٢��0`�ߢq�u0`�ߢq�	�u��٢�-,p�+���S�UW��ÿ����ɿ���g�̿�'�KNHɿR��c�����꟧2�� ��v���?�&6�ld�?Ɨ}B��?���4��?�����?K�Ɓ�?�@��o\�?U�]�_�?�aauw?�aauw?O�]�_�?�@��o\�?K�Ɓ�?�����?���4��?Ɨ}B��?�&6�ld�?��v���?���v�����&6�ldǿ�ŗ}B�п���4�ҿ����ѿG�Ɓʿ�@��o\��C�]�_���aauw��aauw�2�]�_���@��o\��I�Ɓʿ���ѿ���4�ҿƗ}B�п�&6�ldǿ��v�����WMK��?�V��?��_��"�?��t�LS�?h���,��?��e��C�?������?YO���?�o,s�|?|o,s�|?YO���?������?��e��C�?h���,��?��t�LS�?��_��"�?�V��?�WMK��?�WMK����V�̿��_��"Կ��t�LSֿj���,�Կ��e��CϿ�����¿YO�����o,s�|�eo,s�|��XO���������¿��e��CϿg���,�Կ��t�LSֿ��_��"Կ�V�̿�WMK���d�x+t�?�Pl,���?��SF��?�/�]�?3Ib�?�]gK��?�|�d&�?�n�/�b�?v�ck�~?]�ck�~?�n�/�b�?�|�d&�?�]gK��?2Ib�?�/�]�?��SF��?�Pl,���?f�x+t�?T�x+t絿�Pl,��ο ��SF�տ�/�]ؿ5Ibֿ�]gK�ѿ�|�d&Ŀ�n�/�b����ck�~�D�ck�~��n�/�b���|�d&Ŀ�]gK�ѿ1Ibֿ�/�]ؿ��SF�տ�Pl,��οo�x+t絿{o���a?�}��`?�(j:1�Z?nobA~�S?�e� N9I?7�`3��:?nB�o�%?Bu��v�?�O�<�.�>�O�<�.��;u��v��iB�o�%�7�`3��:��e� N9I�qobA~�S��(j:1�Z��}��`�{o���a�{o���a��}��`��(j:1�Z�sobA~�S��e� N9I�&�`3��:��B�o�%�-u��v���O�<�.��hO�<�.�>u��v�?UB�o�%?3�`3��:?�e� N9I?jobA~�S?�(j:1�Z?�}��`?{o���a?�H�׵�?�!�8�o�?٣�U\�?�۸�Ɩ�?u'>��c�?R�{�%�s?\Y�oh_?F(\MC�<?9
��c��>"
��c��8(\MC�<�VY�oh_�R�{�%�s�s'>��c���۸�Ɩ��٣�U\���!�8�o���H�׵���H�׵���!�8�o��ܣ�U\���۸�Ɩ��w'>��c��G�{�%�s�|Y�oh_�'(\MC�<�K
��c��
��c��>(\MC�<?5Y�oh_?O�{�%�s?n'>��c�?�۸�Ɩ�?ң�U\�?�!�8�o�?�H�׵�?wB�l[1�?@䝙l$�?b��v�?��0!ʩ?H[����?`f0+��?_�e U|?�i�Z?�g�`��?ig�`����i�Z��^�e U|�`f0+���E[�������0!ʩ�b��v��@䝙l$��wB�l[1��wB�l[1��B䝙l$��h��v�� �0!ʩ�J[�����Sf0+��� _�e U|��i�Z��g�`���Dg�`��?�i�Z?�^�e U|?\f0+��?C[����?��0!ʩ?^��v�?@䝙l$�?wB�l[1�?FT���? H�&zM�?��cd �?�MD��?~&^gn�?�"�4�Z�?M��H�?� f[	�m?�j(�H"?uj(�H"�� f[	�m�J��H���"�4�Z��{&^gn���MD�����cd Ŀ H�&zMȿFT��ʿFT��ʿ H�&zMȿ��cd Ŀ�MD����&^gn���"�4�Z��]��H��� f[	�m��j(�H"�^j(�H"?� f[	�m?;��H�?�"�4�Z�?w&^gn�?�MD��?��cd �? H�&zM�?FT���?�u6aE��?�br�ع�?�<g@��?�sKce��?1?����?A6%�<2�?#L�{]�?� ���z?A��f�W0?,��f�W0�� ���z�L�{]��A6%�<2��.?������sKce�ʿ�<g@�ѿ�br�عտ�u6aE�׿�u6aE�׿�br�عտ�<g@�ѿ�sKce�ʿ1?�����56%�<2��>L�{]��� ���z�P��f�W0���f�W0?v ���z? L�{]�??6%�<2�?*?����?�sKce��?�<g@��?�br�ع�?�u6aE��?x,#%5 �?"�=��h�?��N���?�E�)�?�����?�J2�|�?*�m��?�|	�4�?`��� �8?D��� �8��|	�4��*�m����J2�|������ɿ��E�)Կ��N��ۿ"�=��h�x,#%5 �x,#%5 �$�=��h���N��ۿ��E�)Կ����ɿ�J2�|��+*�m���t|	�4��|��� �8�$��� �8?^|	�4�?�)�m��?�J2�|�?�����?�E�)�?��N���?"�=��h�?x,#%5 �?o��*��?�N��F��?�x0�A�?̪{fj��?�����?F?��G�?�"؂q?�?��uB1ߊ?|�k��j@?f�k��j@�}�uB1ߊ��"؂q?��F?��G¿���ѿϪ{fj�ڿ�x0�A⿾N��F��o��*��o��*���N��F�志x0�A�Ҫ{fj�ڿ����ѿ<?��G¿#؂q?��j�uB1ߊ���k��j@�S�k��j@?M�uB1ߊ?�"؂q?�?B?��G�?}����?Ū{fj��?�x0�A�?�N��F��?o��*��?Ha��m��?a0rr#�?3e��U��?г�jK��?i���,��?����P��?2�S����?PjV�?�g����C?�g����C�
PjV��.�S���������P�ſg���,�ԿԳ�jK�߿3e��U��a0rr#�Ha��m��Ha��m��a0rr#�:e��U��س�jK�߿l���,�Կ{���P�ſD�S������PjV���g����C��g����C?�PjV�?�S����?����P��?d���,��?ȳ�jK��?/e��U��?a0rr#�?Ha��m��?7�Ȯ�K�?�q�[��?��ʐ�?(�-�Zf�?4Ib�?d�[����?5����?�hD�	��?�gY��uE?�gY��uE��hD�	���1�����d�[���ǿ1Ibֿ*�-�Zf���ʐ��q�[��7�Ȯ�K�7�Ȯ�K��q�[����ʐ�,�-�Zf�7IbֿT�[���ǿF�����xhD�	����gY��uE�ggY��uE?ehD�	��?����?`�[����?.Ib�?$�-�Zf�?��ʐ�?�q�[��?7�Ȯ�K�?
�
Const_1Const*
_output_shapes
:	�*
dtype0*�
value�B�	�"�ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?ĀH&��^?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?����'�q?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?l�LE`Q{?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?��!�`
�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?%=�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?l��(�"�?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?9��,@��?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?t>z6\�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?�#��:�?
`
Const_2Const*
_output_shapes

:
*
dtype0*!
valueB
2        

NoOpNoOp
�|
Const_3Const"/device:CPU:0*
_output_shapes
: *
dtype0*�{
value�{B�{ B�{
�

core_model
	optimizer
loss
	variables
trainable_variables
regularization_losses
	keras_api

signatures
�
	layer-0

layer_with_weights-0

layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer-4
layer_with_weights-3
layer-5
layer-6
layer-7
layer_with_weights-4
layer-8
layer-9
layer-10
layer_with_weights-5
layer-11
layer-12
layer-13
layer_with_weights-6
layer-14
layer-15
layer-16
layer_with_weights-7
layer-17
layer-18
layer-19
layer_with_weights-8
layer-20
layer-21
layer_with_weights-9
layer-22
 	variables
!trainable_variables
"regularization_losses
#	keras_api
�
$iter

%beta_1

&beta_2
	'decay
(learning_rate+m�,m�-m�.m�/m�0m�1m�2m�3m�4m�5m�6m�7m�8m�9m�:m�+v�,v�-v�.v�/v�0v�1v�2v�3v�4v�5v�6v�7v�8v�9v�:v�
 
�
)0
*1
+2
,3
-4
.5
/6
07
18
29
310
411
512
613
714
815
916
:17
v
+0
,1
-2
.3
/4
05
16
27
38
49
510
611
712
813
914
:15
 
�
;non_trainable_variables

<layers
=metrics
>layer_regularization_losses
?layer_metrics
	variables
trainable_variables
regularization_losses
 
 
Z
)mu
@	variables
Atrainable_variables
Bregularization_losses
C	keras_api
b
*
ev_cov_mat
D	variables
Etrainable_variables
Fregularization_losses
G	keras_api
h

+kernel
,bias
H	variables
Itrainable_variables
Jregularization_losses
K	keras_api

L	keras_api
h

-kernel
.bias
M	variables
Ntrainable_variables
Oregularization_losses
P	keras_api
R
Q	variables
Rtrainable_variables
Sregularization_losses
T	keras_api

U	keras_api
h

/kernel
0bias
V	variables
Wtrainable_variables
Xregularization_losses
Y	keras_api
R
Z	variables
[trainable_variables
\regularization_losses
]	keras_api

^	keras_api
h

1kernel
2bias
_	variables
`trainable_variables
aregularization_losses
b	keras_api
R
c	variables
dtrainable_variables
eregularization_losses
f	keras_api

g	keras_api
h

3kernel
4bias
h	variables
itrainable_variables
jregularization_losses
k	keras_api
R
l	variables
mtrainable_variables
nregularization_losses
o	keras_api

p	keras_api
h

5kernel
6bias
q	variables
rtrainable_variables
sregularization_losses
t	keras_api
R
u	variables
vtrainable_variables
wregularization_losses
x	keras_api

y	keras_api
h

7kernel
8bias
z	variables
{trainable_variables
|regularization_losses
}	keras_api
T
~	variables
trainable_variables
�regularization_losses
�	keras_api
l

9kernel
:bias
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�
)0
*1
+2
,3
-4
.5
/6
07
18
29
310
411
512
613
714
815
916
:17
v
+0
,1
-2
.3
/4
05
16
27
38
49
510
611
712
813
914
:15
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
 	variables
!trainable_variables
"regularization_losses
HF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
ZX
VARIABLE_VALUEAdam/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE
DB
VARIABLE_VALUEVariable&variables/0/.ATTRIBUTES/VARIABLE_VALUE
FD
VARIABLE_VALUE
Variable_1&variables/1/.ATTRIBUTES/VARIABLE_VALUE
NL
VARIABLE_VALUElayer_input/kernel&variables/2/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUElayer_input/bias&variables/3/.ATTRIBUTES/VARIABLE_VALUE
RP
VARIABLE_VALUEblock_0_layer_0/kernel&variables/4/.ATTRIBUTES/VARIABLE_VALUE
PN
VARIABLE_VALUEblock_0_layer_0/bias&variables/5/.ATTRIBUTES/VARIABLE_VALUE
RP
VARIABLE_VALUEblock_1_layer_0/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE
PN
VARIABLE_VALUEblock_1_layer_0/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE
RP
VARIABLE_VALUEblock_2_layer_0/kernel&variables/8/.ATTRIBUTES/VARIABLE_VALUE
PN
VARIABLE_VALUEblock_2_layer_0/bias&variables/9/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEblock_3_layer_0/kernel'variables/10/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEblock_3_layer_0/bias'variables/11/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEblock_4_layer_0/kernel'variables/12/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEblock_4_layer_0/bias'variables/13/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEblock_5_layer_0/kernel'variables/14/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEblock_5_layer_0/bias'variables/15/.ATTRIBUTES/VARIABLE_VALUE
PN
VARIABLE_VALUEdense_output/kernel'variables/16/.ATTRIBUTES/VARIABLE_VALUE
NL
VARIABLE_VALUEdense_output/bias'variables/17/.ATTRIBUTES/VARIABLE_VALUE

)0
*1

0
P
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
 
 

)0
 
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
@	variables
Atrainable_variables
Bregularization_losses

*0
 
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
D	variables
Etrainable_variables
Fregularization_losses

+0
,1

+0
,1
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
H	variables
Itrainable_variables
Jregularization_losses
 

-0
.1

-0
.1
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
M	variables
Ntrainable_variables
Oregularization_losses
 
 
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
Q	variables
Rtrainable_variables
Sregularization_losses
 

/0
01

/0
01
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
V	variables
Wtrainable_variables
Xregularization_losses
 
 
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
Z	variables
[trainable_variables
\regularization_losses
 

10
21

10
21
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
_	variables
`trainable_variables
aregularization_losses
 
 
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
c	variables
dtrainable_variables
eregularization_losses
 

30
41

30
41
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
h	variables
itrainable_variables
jregularization_losses
 
 
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
l	variables
mtrainable_variables
nregularization_losses
 

50
61

50
61
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
q	variables
rtrainable_variables
sregularization_losses
 
 
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
u	variables
vtrainable_variables
wregularization_losses
 

70
81

70
81
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
z	variables
{trainable_variables
|regularization_losses
 
 
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
~	variables
trainable_variables
�regularization_losses

90
:1

90
:1
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses

)0
*1
�
	0

1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
 
 
 
8

�total

�count
�	variables
�	keras_api
8

�total

�count
�	variables
�	keras_api
8

�total

�count
�	variables
�	keras_api
8

�total

�count
�	variables
�	keras_api
I

�total

�count
�
_fn_kwargs
�	variables
�	keras_api
I

�total

�count
�
_fn_kwargs
�	variables
�	keras_api
I

�total

�count
�
_fn_kwargs
�	variables
�	keras_api
I

�total

�count
�
_fn_kwargs
�	variables
�	keras_api
I

�total

�count
�
_fn_kwargs
�	variables
�	keras_api
I

�total

�count
�
_fn_kwargs
�	variables
�	keras_api

)0
 
 
 
 

*0
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_14keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_14keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_24keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_24keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUE

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_34keras_api/metrics/3/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_34keras_api/metrics/3/count/.ATTRIBUTES/VARIABLE_VALUE

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_44keras_api/metrics/4/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_44keras_api/metrics/4/count/.ATTRIBUTES/VARIABLE_VALUE
 

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_54keras_api/metrics/5/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_54keras_api/metrics/5/count/.ATTRIBUTES/VARIABLE_VALUE
 

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_64keras_api/metrics/6/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_64keras_api/metrics/6/count/.ATTRIBUTES/VARIABLE_VALUE
 

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_74keras_api/metrics/7/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_74keras_api/metrics/7/count/.ATTRIBUTES/VARIABLE_VALUE
 

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_84keras_api/metrics/8/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_84keras_api/metrics/8/count/.ATTRIBUTES/VARIABLE_VALUE
 

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_94keras_api/metrics/9/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_94keras_api/metrics/9/count/.ATTRIBUTES/VARIABLE_VALUE
 

�0
�1

�	variables
qo
VARIABLE_VALUEAdam/layer_input/kernel/mBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/layer_input/bias/mBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEAdam/block_0_layer_0/kernel/mBvariables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/block_0_layer_0/bias/mBvariables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEAdam/block_1_layer_0/kernel/mBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/block_1_layer_0/bias/mBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEAdam/block_2_layer_0/kernel/mBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/block_2_layer_0/bias/mBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/block_3_layer_0/kernel/mCvariables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/block_3_layer_0/bias/mCvariables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/block_4_layer_0/kernel/mCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/block_4_layer_0/bias/mCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/block_5_layer_0/kernel/mCvariables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/block_5_layer_0/bias/mCvariables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/dense_output/kernel/mCvariables/16/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/dense_output/bias/mCvariables/17/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/layer_input/kernel/vBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/layer_input/bias/vBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEAdam/block_0_layer_0/kernel/vBvariables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/block_0_layer_0/bias/vBvariables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEAdam/block_1_layer_0/kernel/vBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/block_1_layer_0/bias/vBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEAdam/block_2_layer_0/kernel/vBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/block_2_layer_0/bias/vBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/block_3_layer_0/kernel/vCvariables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/block_3_layer_0/bias/vCvariables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/block_4_layer_0/kernel/vCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/block_4_layer_0/bias/vCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/block_5_layer_0/kernel/vCvariables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/block_5_layer_0/bias/vCvariables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/dense_output/kernel/vCvariables/16/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/dense_output/bias/vCvariables/17/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
z
serving_default_input_1Placeholder*'
_output_shapes
:���������	*
dtype0*
shape:���������	
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1Variable
Variable_1layer_input/kernellayer_input/biasblock_0_layer_0/kernelblock_0_layer_0/biasblock_1_layer_0/kernelblock_1_layer_0/biasblock_2_layer_0/kernelblock_2_layer_0/biasblock_3_layer_0/kernelblock_3_layer_0/biasblock_4_layer_0/kernelblock_4_layer_0/biasblock_5_layer_0/kernelblock_5_layer_0/biasdense_output/kerneldense_output/biasConstConst_1Const_2*!
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������	:���������	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� */
f*R(
&__inference_signature_wrapper_93978582
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenameAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOpVariable/Read/ReadVariableOpVariable_1/Read/ReadVariableOp&layer_input/kernel/Read/ReadVariableOp$layer_input/bias/Read/ReadVariableOp*block_0_layer_0/kernel/Read/ReadVariableOp(block_0_layer_0/bias/Read/ReadVariableOp*block_1_layer_0/kernel/Read/ReadVariableOp(block_1_layer_0/bias/Read/ReadVariableOp*block_2_layer_0/kernel/Read/ReadVariableOp(block_2_layer_0/bias/Read/ReadVariableOp*block_3_layer_0/kernel/Read/ReadVariableOp(block_3_layer_0/bias/Read/ReadVariableOp*block_4_layer_0/kernel/Read/ReadVariableOp(block_4_layer_0/bias/Read/ReadVariableOp*block_5_layer_0/kernel/Read/ReadVariableOp(block_5_layer_0/bias/Read/ReadVariableOp'dense_output/kernel/Read/ReadVariableOp%dense_output/bias/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal_2/Read/ReadVariableOpcount_2/Read/ReadVariableOptotal_3/Read/ReadVariableOpcount_3/Read/ReadVariableOptotal_4/Read/ReadVariableOpcount_4/Read/ReadVariableOptotal_5/Read/ReadVariableOpcount_5/Read/ReadVariableOptotal_6/Read/ReadVariableOpcount_6/Read/ReadVariableOptotal_7/Read/ReadVariableOpcount_7/Read/ReadVariableOptotal_8/Read/ReadVariableOpcount_8/Read/ReadVariableOptotal_9/Read/ReadVariableOpcount_9/Read/ReadVariableOp-Adam/layer_input/kernel/m/Read/ReadVariableOp+Adam/layer_input/bias/m/Read/ReadVariableOp1Adam/block_0_layer_0/kernel/m/Read/ReadVariableOp/Adam/block_0_layer_0/bias/m/Read/ReadVariableOp1Adam/block_1_layer_0/kernel/m/Read/ReadVariableOp/Adam/block_1_layer_0/bias/m/Read/ReadVariableOp1Adam/block_2_layer_0/kernel/m/Read/ReadVariableOp/Adam/block_2_layer_0/bias/m/Read/ReadVariableOp1Adam/block_3_layer_0/kernel/m/Read/ReadVariableOp/Adam/block_3_layer_0/bias/m/Read/ReadVariableOp1Adam/block_4_layer_0/kernel/m/Read/ReadVariableOp/Adam/block_4_layer_0/bias/m/Read/ReadVariableOp1Adam/block_5_layer_0/kernel/m/Read/ReadVariableOp/Adam/block_5_layer_0/bias/m/Read/ReadVariableOp.Adam/dense_output/kernel/m/Read/ReadVariableOp,Adam/dense_output/bias/m/Read/ReadVariableOp-Adam/layer_input/kernel/v/Read/ReadVariableOp+Adam/layer_input/bias/v/Read/ReadVariableOp1Adam/block_0_layer_0/kernel/v/Read/ReadVariableOp/Adam/block_0_layer_0/bias/v/Read/ReadVariableOp1Adam/block_1_layer_0/kernel/v/Read/ReadVariableOp/Adam/block_1_layer_0/bias/v/Read/ReadVariableOp1Adam/block_2_layer_0/kernel/v/Read/ReadVariableOp/Adam/block_2_layer_0/bias/v/Read/ReadVariableOp1Adam/block_3_layer_0/kernel/v/Read/ReadVariableOp/Adam/block_3_layer_0/bias/v/Read/ReadVariableOp1Adam/block_4_layer_0/kernel/v/Read/ReadVariableOp/Adam/block_4_layer_0/bias/v/Read/ReadVariableOp1Adam/block_5_layer_0/kernel/v/Read/ReadVariableOp/Adam/block_5_layer_0/bias/v/Read/ReadVariableOp.Adam/dense_output/kernel/v/Read/ReadVariableOp,Adam/dense_output/bias/v/Read/ReadVariableOpConst_3*X
TinQ
O2M	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� **
f%R#
!__inference__traced_save_93980800
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filename	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_rateVariable
Variable_1layer_input/kernellayer_input/biasblock_0_layer_0/kernelblock_0_layer_0/biasblock_1_layer_0/kernelblock_1_layer_0/biasblock_2_layer_0/kernelblock_2_layer_0/biasblock_3_layer_0/kernelblock_3_layer_0/biasblock_4_layer_0/kernelblock_4_layer_0/biasblock_5_layer_0/kernelblock_5_layer_0/biasdense_output/kerneldense_output/biastotalcounttotal_1count_1total_2count_2total_3count_3total_4count_4total_5count_5total_6count_6total_7count_7total_8count_8total_9count_9Adam/layer_input/kernel/mAdam/layer_input/bias/mAdam/block_0_layer_0/kernel/mAdam/block_0_layer_0/bias/mAdam/block_1_layer_0/kernel/mAdam/block_1_layer_0/bias/mAdam/block_2_layer_0/kernel/mAdam/block_2_layer_0/bias/mAdam/block_3_layer_0/kernel/mAdam/block_3_layer_0/bias/mAdam/block_4_layer_0/kernel/mAdam/block_4_layer_0/bias/mAdam/block_5_layer_0/kernel/mAdam/block_5_layer_0/bias/mAdam/dense_output/kernel/mAdam/dense_output/bias/mAdam/layer_input/kernel/vAdam/layer_input/bias/vAdam/block_0_layer_0/kernel/vAdam/block_0_layer_0/bias/vAdam/block_1_layer_0/kernel/vAdam/block_1_layer_0/bias/vAdam/block_2_layer_0/kernel/vAdam/block_2_layer_0/bias/vAdam/block_3_layer_0/kernel/vAdam/block_3_layer_0/bias/vAdam/block_4_layer_0/kernel/vAdam/block_4_layer_0/bias/vAdam/block_5_layer_0/kernel/vAdam/block_5_layer_0/bias/vAdam/dense_output/kernel/vAdam/dense_output/bias/v*W
TinP
N2L*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *-
f(R&
$__inference__traced_restore_93981035֏5
��
�
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93975642

inputs'
mean_shift_layer_93975265:	.
decorrelation_layer_93975276:		'
layer_input_93975302:		�#
layer_input_93975304:	�,
block_0_layer_0_93975331:
��'
block_0_layer_0_93975333:	�,
block_1_layer_0_93975368:
��'
block_1_layer_0_93975370:	�,
block_2_layer_0_93975405:
��'
block_2_layer_0_93975407:	�,
block_3_layer_0_93975442:
��'
block_3_layer_0_93975444:	�,
block_4_layer_0_93975479:
��'
block_4_layer_0_93975481:	�,
block_5_layer_0_93975516:
��'
block_5_layer_0_93975518:	�(
dense_output_93975546:	�#
dense_output_93975548:
identity��'block_0_layer_0/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_1_layer_0/StatefulPartitionedCall�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_2_layer_0/StatefulPartitionedCall�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_3_layer_0/StatefulPartitionedCall�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_4_layer_0/StatefulPartitionedCall�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_5_layer_0/StatefulPartitionedCall�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�+decorrelation_layer/StatefulPartitionedCall�$dense_output/StatefulPartitionedCall�5dense_output/kernel/Regularizer/Square/ReadVariableOp�#layer_input/StatefulPartitionedCall�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�(mean_shift_layer/StatefulPartitionedCall�
(mean_shift_layer/StatefulPartitionedCallStatefulPartitionedCallinputsmean_shift_layer_93975265*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *W
fRRP
N__inference_mean_shift_layer_layer_call_and_return_conditional_losses_93975264�
+decorrelation_layer/StatefulPartitionedCallStatefulPartitionedCall1mean_shift_layer/StatefulPartitionedCall:output:0decorrelation_layer_93975276*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *Z
fURS
Q__inference_decorrelation_layer_layer_call_and_return_conditional_losses_93975275�
#layer_input/StatefulPartitionedCallStatefulPartitionedCall4decorrelation_layer/StatefulPartitionedCall:output:0layer_input_93975302layer_input_93975304*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *R
fMRK
I__inference_layer_input_layer_call_and_return_conditional_losses_93975301�
tf.math.softplus/SoftplusSoftplus,layer_input/StatefulPartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_0_layer_0/StatefulPartitionedCallStatefulPartitionedCall'tf.math.softplus/Softplus:activations:0block_0_layer_0_93975331block_0_layer_0_93975333*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_0_layer_0_layer_call_and_return_conditional_losses_93975330�
add/PartitionedCallPartitionedCall,layer_input/StatefulPartitionedCall:output:00block_0_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__inference_add_layer_call_and_return_conditional_losses_93975342x
tf.math.softplus_1/SoftplusSoftplusadd/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_1_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_1/Softplus:activations:0block_1_layer_0_93975368block_1_layer_0_93975370*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_1_layer_0_layer_call_and_return_conditional_losses_93975367�
add_1/PartitionedCallPartitionedCalladd/PartitionedCall:output:00block_1_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_1_layer_call_and_return_conditional_losses_93975379z
tf.math.softplus_2/SoftplusSoftplusadd_1/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_2_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_2/Softplus:activations:0block_2_layer_0_93975405block_2_layer_0_93975407*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_2_layer_0_layer_call_and_return_conditional_losses_93975404�
add_2/PartitionedCallPartitionedCalladd_1/PartitionedCall:output:00block_2_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_2_layer_call_and_return_conditional_losses_93975416z
tf.math.softplus_3/SoftplusSoftplusadd_2/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_3_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_3/Softplus:activations:0block_3_layer_0_93975442block_3_layer_0_93975444*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_3_layer_0_layer_call_and_return_conditional_losses_93975441�
add_3/PartitionedCallPartitionedCalladd_2/PartitionedCall:output:00block_3_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_3_layer_call_and_return_conditional_losses_93975453z
tf.math.softplus_4/SoftplusSoftplusadd_3/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_4_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_4/Softplus:activations:0block_4_layer_0_93975479block_4_layer_0_93975481*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_4_layer_0_layer_call_and_return_conditional_losses_93975478�
add_4/PartitionedCallPartitionedCalladd_3/PartitionedCall:output:00block_4_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_4_layer_call_and_return_conditional_losses_93975490z
tf.math.softplus_5/SoftplusSoftplusadd_4/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_5_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_5/Softplus:activations:0block_5_layer_0_93975516block_5_layer_0_93975518*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_5_layer_0_layer_call_and_return_conditional_losses_93975515�
add_5/PartitionedCallPartitionedCalladd_4/PartitionedCall:output:00block_5_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_5_layer_call_and_return_conditional_losses_93975527�
$dense_output/StatefulPartitionedCallStatefulPartitionedCalladd_5/PartitionedCall:output:0dense_output_93975546dense_output_93975548*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *S
fNRL
J__inference_dense_output_layer_call_and_return_conditional_losses_93975545�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93975302*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93975304*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93975331* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93975333*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93975368* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93975370*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93975405* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93975407*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93975442* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93975444*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93975479* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93975481*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93975516* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93975518*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_output_93975546*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: |
IdentityIdentity-dense_output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������

NoOpNoOp(^block_0_layer_0/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_1_layer_0/StatefulPartitionedCall7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_2_layer_0/StatefulPartitionedCall7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_3_layer_0/StatefulPartitionedCall7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_4_layer_0/StatefulPartitionedCall7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_5_layer_0/StatefulPartitionedCall7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp,^decorrelation_layer/StatefulPartitionedCall%^dense_output/StatefulPartitionedCall6^dense_output/kernel/Regularizer/Square/ReadVariableOp$^layer_input/StatefulPartitionedCall3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp)^mean_shift_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 2R
'block_0_layer_0/StatefulPartitionedCall'block_0_layer_0/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_1_layer_0/StatefulPartitionedCall'block_1_layer_0/StatefulPartitionedCall2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_2_layer_0/StatefulPartitionedCall'block_2_layer_0/StatefulPartitionedCall2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_3_layer_0/StatefulPartitionedCall'block_3_layer_0/StatefulPartitionedCall2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_4_layer_0/StatefulPartitionedCall'block_4_layer_0/StatefulPartitionedCall2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_5_layer_0/StatefulPartitionedCall'block_5_layer_0/StatefulPartitionedCall2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2Z
+decorrelation_layer/StatefulPartitionedCall+decorrelation_layer/StatefulPartitionedCall2L
$dense_output/StatefulPartitionedCall$dense_output/StatefulPartitionedCall2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2J
#layer_input/StatefulPartitionedCall#layer_input/StatefulPartitionedCall2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp2T
(mean_shift_layer/StatefulPartitionedCall(mean_shift_layer/StatefulPartitionedCall:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
K__forward_block_2_layer_0_layer_call_and_return_conditional_losses_93976653
inputs_02
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity
matmul_readvariableop

inputs��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0l
MatMulMatMulinputs_0MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"6
matmul_readvariableopMatMul/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : *}
backward_function_nameca__inference___backward_block_2_layer_0_layer_call_and_return_conditional_losses_93976641_9397665420
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
o
C__inference_add_5_layer_call_and_return_conditional_losses_93980351
inputs_0
inputs_1
identityS
addAddV2inputs_0inputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
��
�
h__inference___backward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976445_93976848
placeholder\
Xgradients_dense_output_statefulpartitionedcall_grad_dense_output_statefulpartitionedcall^
Zgradients_dense_output_statefulpartitionedcall_grad_dense_output_statefulpartitionedcall_1>
:gradients_add_5_partitionedcall_grad_add_5_partitionedcall@
<gradients_add_5_partitionedcall_grad_add_5_partitionedcall_1b
^gradients_block_5_layer_0_statefulpartitionedcall_grad_block_5_layer_0_statefulpartitionedcalld
`gradients_block_5_layer_0_statefulpartitionedcall_grad_block_5_layer_0_statefulpartitionedcall_1L
Hgradients_tf_math_softplus_5_softplus_grad_sigmoid_add_4_partitionedcall>
:gradients_add_4_partitionedcall_grad_add_4_partitionedcall@
<gradients_add_4_partitionedcall_grad_add_4_partitionedcall_1b
^gradients_block_4_layer_0_statefulpartitionedcall_grad_block_4_layer_0_statefulpartitionedcalld
`gradients_block_4_layer_0_statefulpartitionedcall_grad_block_4_layer_0_statefulpartitionedcall_1L
Hgradients_tf_math_softplus_4_softplus_grad_sigmoid_add_3_partitionedcall>
:gradients_add_3_partitionedcall_grad_add_3_partitionedcall@
<gradients_add_3_partitionedcall_grad_add_3_partitionedcall_1b
^gradients_block_3_layer_0_statefulpartitionedcall_grad_block_3_layer_0_statefulpartitionedcalld
`gradients_block_3_layer_0_statefulpartitionedcall_grad_block_3_layer_0_statefulpartitionedcall_1L
Hgradients_tf_math_softplus_3_softplus_grad_sigmoid_add_2_partitionedcall>
:gradients_add_2_partitionedcall_grad_add_2_partitionedcall@
<gradients_add_2_partitionedcall_grad_add_2_partitionedcall_1b
^gradients_block_2_layer_0_statefulpartitionedcall_grad_block_2_layer_0_statefulpartitionedcalld
`gradients_block_2_layer_0_statefulpartitionedcall_grad_block_2_layer_0_statefulpartitionedcall_1L
Hgradients_tf_math_softplus_2_softplus_grad_sigmoid_add_1_partitionedcall>
:gradients_add_1_partitionedcall_grad_add_1_partitionedcall@
<gradients_add_1_partitionedcall_grad_add_1_partitionedcall_1b
^gradients_block_1_layer_0_statefulpartitionedcall_grad_block_1_layer_0_statefulpartitionedcalld
`gradients_block_1_layer_0_statefulpartitionedcall_grad_block_1_layer_0_statefulpartitionedcall_1J
Fgradients_tf_math_softplus_1_softplus_grad_sigmoid_add_partitionedcall:
6gradients_add_partitionedcall_grad_add_partitionedcall<
8gradients_add_partitionedcall_grad_add_partitionedcall_1b
^gradients_block_0_layer_0_statefulpartitionedcall_grad_block_0_layer_0_statefulpartitionedcalld
`gradients_block_0_layer_0_statefulpartitionedcall_grad_block_0_layer_0_statefulpartitionedcall_1X
Tgradients_tf_math_softplus_softplus_grad_sigmoid_layer_input_statefulpartitionedcallZ
Vgradients_layer_input_statefulpartitionedcall_grad_layer_input_statefulpartitionedcall\
Xgradients_layer_input_statefulpartitionedcall_grad_layer_input_statefulpartitionedcall_1j
fgradients_decorrelation_layer_statefulpartitionedcall_grad_decorrelation_layer_statefulpartitionedcalll
hgradients_decorrelation_layer_statefulpartitionedcall_grad_decorrelation_layer_statefulpartitionedcall_1d
`gradients_mean_shift_layer_statefulpartitionedcall_grad_mean_shift_layer_statefulpartitionedcallf
bgradients_mean_shift_layer_statefulpartitionedcall_grad_mean_shift_layer_statefulpartitionedcall_1
identity

identity_1

identity_2

identity_3

identity_4

identity_5

identity_6

identity_7

identity_8

identity_9
identity_10
identity_11
identity_12
identity_13
identity_14
identity_15
identity_16
identity_17
identity_18^
gradients/grad_ys_0Identityplaceholder*
T0*'
_output_shapes
:����������
Cgradients/dense_output/StatefulPartitionedCall_grad/PartitionedCallPartitionedCallgradients/grad_ys_0:output:0Xgradients_dense_output_statefulpartitionedcall_grad_dense_output_statefulpartitionedcallZgradients_dense_output_statefulpartitionedcall_grad_dense_output_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *9
_output_shapes'
%:����������:	�:* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *g
fbR`
^__inference___backward_dense_output_layer_call_and_return_conditional_losses_93976449_93976462�
4gradients/add_5/PartitionedCall_grad/PartitionedCallPartitionedCallLgradients/dense_output/StatefulPartitionedCall_grad/PartitionedCall:output:0:gradients_add_5_partitionedcall_grad_add_5_partitionedcall<gradients_add_5_partitionedcall_grad_add_5_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_5_layer_call_and_return_conditional_losses_93976470_93976487�
Fgradients/block_5_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_5/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_5_layer_0_statefulpartitionedcall_grad_block_5_layer_0_statefulpartitionedcall`gradients_block_5_layer_0_statefulpartitionedcall_grad_block_5_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_5_layer_0_layer_call_and_return_conditional_losses_93976494_93976507�
2gradients/tf.math.softplus_5/Softplus_grad/SigmoidSigmoidHgradients_tf_math_softplus_5_softplus_grad_sigmoid_add_4_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_5/Softplus_grad/mulMulOgradients/block_5_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_5/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddNAddN=gradients/add_5/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_5/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_5/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
4gradients/add_4/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN:sum:0:gradients_add_4_partitionedcall_grad_add_4_partitionedcall<gradients_add_4_partitionedcall_grad_add_4_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_4_layer_call_and_return_conditional_losses_93976519_93976536�
Fgradients/block_4_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_4/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_4_layer_0_statefulpartitionedcall_grad_block_4_layer_0_statefulpartitionedcall`gradients_block_4_layer_0_statefulpartitionedcall_grad_block_4_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_4_layer_0_layer_call_and_return_conditional_losses_93976543_93976556�
2gradients/tf.math.softplus_4/Softplus_grad/SigmoidSigmoidHgradients_tf_math_softplus_4_softplus_grad_sigmoid_add_3_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_4/Softplus_grad/mulMulOgradients/block_4_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_4/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_1AddN=gradients/add_4/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_4/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_4/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
4gradients/add_3/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_1:sum:0:gradients_add_3_partitionedcall_grad_add_3_partitionedcall<gradients_add_3_partitionedcall_grad_add_3_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_3_layer_call_and_return_conditional_losses_93976568_93976585�
Fgradients/block_3_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_3/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_3_layer_0_statefulpartitionedcall_grad_block_3_layer_0_statefulpartitionedcall`gradients_block_3_layer_0_statefulpartitionedcall_grad_block_3_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_3_layer_0_layer_call_and_return_conditional_losses_93976592_93976605�
2gradients/tf.math.softplus_3/Softplus_grad/SigmoidSigmoidHgradients_tf_math_softplus_3_softplus_grad_sigmoid_add_2_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_3/Softplus_grad/mulMulOgradients/block_3_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_3/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_2AddN=gradients/add_3/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_3/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_3/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
4gradients/add_2/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_2:sum:0:gradients_add_2_partitionedcall_grad_add_2_partitionedcall<gradients_add_2_partitionedcall_grad_add_2_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_2_layer_call_and_return_conditional_losses_93976617_93976634�
Fgradients/block_2_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_2/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_2_layer_0_statefulpartitionedcall_grad_block_2_layer_0_statefulpartitionedcall`gradients_block_2_layer_0_statefulpartitionedcall_grad_block_2_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_2_layer_0_layer_call_and_return_conditional_losses_93976641_93976654�
2gradients/tf.math.softplus_2/Softplus_grad/SigmoidSigmoidHgradients_tf_math_softplus_2_softplus_grad_sigmoid_add_1_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_2/Softplus_grad/mulMulOgradients/block_2_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_2/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_3AddN=gradients/add_2/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_2/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_2/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
4gradients/add_1/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_3:sum:0:gradients_add_1_partitionedcall_grad_add_1_partitionedcall<gradients_add_1_partitionedcall_grad_add_1_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_1_layer_call_and_return_conditional_losses_93976666_93976683�
Fgradients/block_1_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_1/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_1_layer_0_statefulpartitionedcall_grad_block_1_layer_0_statefulpartitionedcall`gradients_block_1_layer_0_statefulpartitionedcall_grad_block_1_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_1_layer_0_layer_call_and_return_conditional_losses_93976690_93976703�
2gradients/tf.math.softplus_1/Softplus_grad/SigmoidSigmoidFgradients_tf_math_softplus_1_softplus_grad_sigmoid_add_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_1/Softplus_grad/mulMulOgradients/block_1_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_1/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_4AddN=gradients/add_1/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_1/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_1/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
2gradients/add/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_4:sum:06gradients_add_partitionedcall_grad_add_partitionedcall8gradients_add_partitionedcall_grad_add_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *^
fYRW
U__inference___backward_add_layer_call_and_return_conditional_losses_93976715_93976732�
Fgradients/block_0_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall;gradients/add/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_0_layer_0_statefulpartitionedcall_grad_block_0_layer_0_statefulpartitionedcall`gradients_block_0_layer_0_statefulpartitionedcall_grad_block_0_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_0_layer_0_layer_call_and_return_conditional_losses_93976739_93976752�
0gradients/tf.math.softplus/Softplus_grad/SigmoidSigmoidTgradients_tf_math_softplus_softplus_grad_sigmoid_layer_input_statefulpartitionedcall*
T0*(
_output_shapes
:�����������
,gradients/tf.math.softplus/Softplus_grad/mulMulOgradients/block_0_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:04gradients/tf.math.softplus/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_5AddN;gradients/add/PartitionedCall_grad/PartitionedCall:output:00gradients/tf.math.softplus/Softplus_grad/mul:z:0*
N*
T0*E
_class;
97loc:@gradients/add/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
Bgradients/layer_input/StatefulPartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_5:sum:0Vgradients_layer_input_statefulpartitionedcall_grad_layer_input_statefulpartitionedcallXgradients_layer_input_statefulpartitionedcall_grad_layer_input_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *9
_output_shapes'
%:���������	:		�:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *f
faR_
]__inference___backward_layer_input_layer_call_and_return_conditional_losses_93976764_93976777�
Jgradients/decorrelation_layer/StatefulPartitionedCall_grad/PartitionedCallPartitionedCallKgradients/layer_input/StatefulPartitionedCall_grad/PartitionedCall:output:0fgradients_decorrelation_layer_statefulpartitionedcall_grad_decorrelation_layer_statefulpartitionedcallhgradients_decorrelation_layer_statefulpartitionedcall_grad_decorrelation_layer_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *1
_output_shapes
:���������	:		* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *n
fiRg
e__inference___backward_decorrelation_layer_layer_call_and_return_conditional_losses_93976785_93976796�
Ggradients/mean_shift_layer/StatefulPartitionedCall_grad/PartitionedCallPartitionedCallSgradients/decorrelation_layer/StatefulPartitionedCall_grad/PartitionedCall:output:0`gradients_mean_shift_layer_statefulpartitionedcall_grad_mean_shift_layer_statefulpartitionedcallbgradients_mean_shift_layer_statefulpartitionedcall_grad_mean_shift_layer_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *-
_output_shapes
:���������	:	* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *k
ffRd
b__inference___backward_mean_shift_layer_layer_call_and_return_conditional_losses_93976803_93976821�
IdentityIdentityPgradients/mean_shift_layer/StatefulPartitionedCall_grad/PartitionedCall:output:0*
T0*'
_output_shapes
:���������	�

Identity_1IdentityPgradients/mean_shift_layer/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0*
_output_shapes
:	�

Identity_2IdentitySgradients/decorrelation_layer/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0*
_output_shapes

:		�

Identity_3IdentityKgradients/layer_input/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0*
_output_shapes
:		��

Identity_4IdentityKgradients/layer_input/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��

Identity_5IdentityOgradients/block_0_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���

Identity_6IdentityOgradients/block_0_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��

Identity_7IdentityOgradients/block_1_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���

Identity_8IdentityOgradients/block_1_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��

Identity_9IdentityOgradients/block_2_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���
Identity_10IdentityOgradients/block_2_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��
Identity_11IdentityOgradients/block_3_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���
Identity_12IdentityOgradients/block_3_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��
Identity_13IdentityOgradients/block_4_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���
Identity_14IdentityOgradients/block_4_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��
Identity_15IdentityOgradients/block_5_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���
Identity_16IdentityOgradients/block_5_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��
Identity_17IdentityLgradients/dense_output/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0*
_output_shapes
:	��
Identity_18IdentityLgradients/dense_output/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes
:"
identityIdentity:output:0"!

identity_1Identity_1:output:0"#
identity_10Identity_10:output:0"#
identity_11Identity_11:output:0"#
identity_12Identity_12:output:0"#
identity_13Identity_13:output:0"#
identity_14Identity_14:output:0"#
identity_15Identity_15:output:0"#
identity_16Identity_16:output:0"#
identity_17Identity_17:output:0"#
identity_18Identity_18:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0"!

identity_4Identity_4:output:0"!

identity_5Identity_5:output:0"!

identity_6Identity_6:output:0"!

identity_7Identity_7:output:0"!

identity_8Identity_8:output:0"!

identity_9Identity_9:output:0*(
_construction_contextkEagerRuntime*�
_input_shapes�
�:���������:	�:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:		�:���������	:		:���������	:���������	:	*m
forward_function_nameTR__forward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976847:- )
'
_output_shapes
:���������:%!

_output_shapes
:	�:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.	*
(
_output_shapes
:����������:&
"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:. *
(
_output_shapes
:����������:%!!

_output_shapes
:		�:-")
'
_output_shapes
:���������	:$# 

_output_shapes

:		:-$)
'
_output_shapes
:���������	:-%)
'
_output_shapes
:���������	: &

_output_shapes
:	
�
�
2__inference_block_3_layer_0_layer_call_fn_93980207

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_3_layer_0_layer_call_and_return_conditional_losses_93975441p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93978115
input_1-
resnet_entropy_closure_93977800:	1
resnet_entropy_closure_93977802:		2
resnet_entropy_closure_93977804:		�.
resnet_entropy_closure_93977806:	�3
resnet_entropy_closure_93977808:
��.
resnet_entropy_closure_93977810:	�3
resnet_entropy_closure_93977812:
��.
resnet_entropy_closure_93977814:	�3
resnet_entropy_closure_93977816:
��.
resnet_entropy_closure_93977818:	�3
resnet_entropy_closure_93977820:
��.
resnet_entropy_closure_93977822:	�3
resnet_entropy_closure_93977824:
��.
resnet_entropy_closure_93977826:	�3
resnet_entropy_closure_93977828:
��.
resnet_entropy_closure_93977830:	�2
resnet_entropy_closure_93977832:	�-
resnet_entropy_closure_93977834:
unknown
tensordot_1_b
mul_1_x
identity

identity_1

identity_2��.ResNet_entropy_closure/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�checked�	checked_1�5dense_output/kernel/Regularizer/Square/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�
.ResNet_entropy_closure/StatefulPartitionedCallStatefulPartitionedCallinput_1resnet_entropy_closure_93977800resnet_entropy_closure_93977802resnet_entropy_closure_93977804resnet_entropy_closure_93977806resnet_entropy_closure_93977808resnet_entropy_closure_93977810resnet_entropy_closure_93977812resnet_entropy_closure_93977814resnet_entropy_closure_93977816resnet_entropy_closure_93977818resnet_entropy_closure_93977820resnet_entropy_closure_93977822resnet_entropy_closure_93977824resnet_entropy_closure_93977826resnet_entropy_closure_93977828resnet_entropy_closure_93977830resnet_entropy_closure_93977832resnet_entropy_closure_93977834*
Tin
2*3
Tout+
)2'*
_collective_manager_ids
 *�
_output_shapes�
�:���������:	�:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:		�:���������	:		:���������	:���������	:	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *[
fVRT
R__forward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976847v
ones_like/ShapeShape7ResNet_entropy_closure/StatefulPartitionedCall:output:0*
T0*
_output_shapes
:T
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?w
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:����������
PartitionedCallPartitionedCallones_like:output:07ResNet_entropy_closure/StatefulPartitionedCall:output:17ResNet_entropy_closure/StatefulPartitionedCall:output:27ResNet_entropy_closure/StatefulPartitionedCall:output:37ResNet_entropy_closure/StatefulPartitionedCall:output:47ResNet_entropy_closure/StatefulPartitionedCall:output:57ResNet_entropy_closure/StatefulPartitionedCall:output:67ResNet_entropy_closure/StatefulPartitionedCall:output:77ResNet_entropy_closure/StatefulPartitionedCall:output:87ResNet_entropy_closure/StatefulPartitionedCall:output:98ResNet_entropy_closure/StatefulPartitionedCall:output:108ResNet_entropy_closure/StatefulPartitionedCall:output:118ResNet_entropy_closure/StatefulPartitionedCall:output:128ResNet_entropy_closure/StatefulPartitionedCall:output:138ResNet_entropy_closure/StatefulPartitionedCall:output:148ResNet_entropy_closure/StatefulPartitionedCall:output:158ResNet_entropy_closure/StatefulPartitionedCall:output:168ResNet_entropy_closure/StatefulPartitionedCall:output:178ResNet_entropy_closure/StatefulPartitionedCall:output:188ResNet_entropy_closure/StatefulPartitionedCall:output:198ResNet_entropy_closure/StatefulPartitionedCall:output:208ResNet_entropy_closure/StatefulPartitionedCall:output:218ResNet_entropy_closure/StatefulPartitionedCall:output:228ResNet_entropy_closure/StatefulPartitionedCall:output:238ResNet_entropy_closure/StatefulPartitionedCall:output:248ResNet_entropy_closure/StatefulPartitionedCall:output:258ResNet_entropy_closure/StatefulPartitionedCall:output:268ResNet_entropy_closure/StatefulPartitionedCall:output:278ResNet_entropy_closure/StatefulPartitionedCall:output:288ResNet_entropy_closure/StatefulPartitionedCall:output:298ResNet_entropy_closure/StatefulPartitionedCall:output:308ResNet_entropy_closure/StatefulPartitionedCall:output:318ResNet_entropy_closure/StatefulPartitionedCall:output:328ResNet_entropy_closure/StatefulPartitionedCall:output:338ResNet_entropy_closure/StatefulPartitionedCall:output:348ResNet_entropy_closure/StatefulPartitionedCall:output:358ResNet_entropy_closure/StatefulPartitionedCall:output:368ResNet_entropy_closure/StatefulPartitionedCall:output:378ResNet_entropy_closure/StatefulPartitionedCall:output:38*2
Tin+
)2'*
Tout
2*
_collective_manager_ids
 *�
_output_shapes�
�:���������	:	:		:		�:�:
��:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�:* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *q
flRj
h__inference___backward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976445_93976848g
CastCastPartitionedCall:output:0*

DstT0*

SrcT0*'
_output_shapes
:���������	�
checkedCheckNumericsCast:y:0*
T0*'
_output_shapes
:���������	*d
messageYWinput tensor checking error at alpha = Tensor("Cast:0", shape=(None, 9), dtype=float64)d
checkedandclipped/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped/MinimumMinimumchecked:output:0$checkedandclipped/Minimum/y:output:0*
T0*'
_output_shapes
:���������	\
checkedandclipped/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclippedMaximumcheckedandclipped/Minimum:z:0checkedandclipped/y:output:0*
T0*'
_output_shapes
:���������	d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSliceunknownstrided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:		�*

begin_mask*
end_maskX
Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:X
Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB: T
Tensordot/ShapeShapecheckedandclipped:z:0*
T0*
_output_shapes
:Y
Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2GatherV2Tensordot/Shape:output:0Tensordot/free:output:0 Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2_1GatherV2Tensordot/Shape:output:0Tensordot/axes:output:0"Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:Y
Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: n
Tensordot/ProdProdTensordot/GatherV2:output:0Tensordot/Const:output:0*
T0*
_output_shapes
: [
Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: t
Tensordot/Prod_1ProdTensordot/GatherV2_1:output:0Tensordot/Const_1:output:0*
T0*
_output_shapes
: W
Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concatConcatV2Tensordot/free:output:0Tensordot/axes:output:0Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:y
Tensordot/stackPackTensordot/Prod:output:0Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot/transpose	Transposecheckedandclipped:z:0Tensordot/concat:output:0*
T0*'
_output_shapes
:���������	�
Tensordot/ReshapeReshapeTensordot/transpose:y:0Tensordot/stack:output:0*
T0*0
_output_shapes
:�������������������
Tensordot/MatMulMatMulTensordot/Reshape:output:0strided_slice:output:0*
T0*(
_output_shapes
:����������\
Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�Y
Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concat_1ConcatV2Tensordot/GatherV2:output:0Tensordot/Const_2:output:0 Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
	TensordotReshapeTensordot/MatMul:product:0Tensordot/concat_1:output:0*
T0*(
_output_shapes
:����������Q
ExpExpTensordot:output:0*
T0*(
_output_shapes
:����������Z
Tensordot_1/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_1/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_1/ShapeShapeExp:y:0*
T0*
_output_shapes
:[
Tensordot_1/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2GatherV2Tensordot_1/Shape:output:0Tensordot_1/free:output:0"Tensordot_1/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_1/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2_1GatherV2Tensordot_1/Shape:output:0Tensordot_1/axes:output:0$Tensordot_1/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_1/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_1/ProdProdTensordot_1/GatherV2:output:0Tensordot_1/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_1/Prod_1ProdTensordot_1/GatherV2_1:output:0Tensordot_1/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concatConcatV2Tensordot_1/free:output:0Tensordot_1/axes:output:0 Tensordot_1/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_1/stackPackTensordot_1/Prod:output:0Tensordot_1/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_1/transpose	TransposeExp:y:0Tensordot_1/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_1/ReshapeReshapeTensordot_1/transpose:y:0Tensordot_1/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_1/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_1/transpose_1	Transposetensordot_1_b%Tensordot_1/transpose_1/perm:output:0*
T0*
_output_shapes
:	��
Tensordot_1/MatMulMatMulTensordot_1/Reshape:output:0Tensordot_1/transpose_1:y:0*
T0*'
_output_shapes
:���������]
Tensordot_1/Const_2Const*
_output_shapes
:*
dtype0*
valueB:[
Tensordot_1/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concat_1ConcatV2Tensordot_1/GatherV2:output:0Tensordot_1/Const_2:output:0"Tensordot_1/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_1ReshapeTensordot_1/MatMul:product:0Tensordot_1/concat_1:output:0*
T0*'
_output_shapes
:���������R
LogLogTensordot_1:output:0*
T0*'
_output_shapes
:���������E
NegNegLog:y:0*
T0*'
_output_shapes
:���������M
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :v
concatConcatV2Neg:y:0Cast:y:0concat/axis:output:0*
N*
T0*'
_output_shapes
:���������
�
	checked_1CheckNumericsconcat:output:0^checked*
T0*'
_output_shapes
:���������
*(
messageinput tensor checking errorf
checkedandclipped_1/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped_1/MinimumMinimumchecked_1:output:0&checkedandclipped_1/Minimum/y:output:0*
T0*'
_output_shapes
:���������
^
checkedandclipped_1/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclipped_1Maximumcheckedandclipped_1/Minimum:z:0checkedandclipped_1/y:output:0*
T0*'
_output_shapes
:���������
Z
Tensordot_2/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_2/freeConst*
_output_shapes
:*
dtype0*
valueB: X
Tensordot_2/ShapeShapecheckedandclipped_1:z:0*
T0*
_output_shapes
:[
Tensordot_2/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2GatherV2Tensordot_2/Shape:output:0Tensordot_2/free:output:0"Tensordot_2/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_2/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2_1GatherV2Tensordot_2/Shape:output:0Tensordot_2/axes:output:0$Tensordot_2/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_2/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_2/ProdProdTensordot_2/GatherV2:output:0Tensordot_2/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_2/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_2/Prod_1ProdTensordot_2/GatherV2_1:output:0Tensordot_2/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_2/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concatConcatV2Tensordot_2/free:output:0Tensordot_2/axes:output:0 Tensordot_2/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_2/stackPackTensordot_2/Prod:output:0Tensordot_2/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2/transpose	Transposecheckedandclipped_1:z:0Tensordot_2/concat:output:0*
T0*'
_output_shapes
:���������
�
Tensordot_2/ReshapeReshapeTensordot_2/transpose:y:0Tensordot_2/stack:output:0*
T0*0
_output_shapes
:������������������v
Tensordot_2/MatMulMatMulTensordot_2/Reshape:output:0unknown*
T0*(
_output_shapes
:����������^
Tensordot_2/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�[
Tensordot_2/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concat_1ConcatV2Tensordot_2/GatherV2:output:0Tensordot_2/Const_2:output:0"Tensordot_2/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2ReshapeTensordot_2/MatMul:product:0Tensordot_2/concat_1:output:0*
T0*(
_output_shapes
:����������U
Exp_1ExpTensordot_2:output:0*
T0*(
_output_shapes
:����������W
MulMul	Exp_1:y:0tensordot_1_b*
T0*(
_output_shapes
:����������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSliceunknownstrided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
:	
�*

begin_mask*
end_maskZ
Tensordot_3/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_3/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_3/ShapeShapeMul:z:0*
T0*
_output_shapes
:[
Tensordot_3/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2GatherV2Tensordot_3/Shape:output:0Tensordot_3/free:output:0"Tensordot_3/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_3/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2_1GatherV2Tensordot_3/Shape:output:0Tensordot_3/axes:output:0$Tensordot_3/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_3/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_3/ProdProdTensordot_3/GatherV2:output:0Tensordot_3/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_3/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_3/Prod_1ProdTensordot_3/GatherV2_1:output:0Tensordot_3/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_3/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concatConcatV2Tensordot_3/free:output:0Tensordot_3/axes:output:0 Tensordot_3/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_3/stackPackTensordot_3/Prod:output:0Tensordot_3/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_3/transpose	TransposeMul:z:0Tensordot_3/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_3/ReshapeReshapeTensordot_3/transpose:y:0Tensordot_3/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_3/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_3/transpose_1	Transposestrided_slice_1:output:0%Tensordot_3/transpose_1/perm:output:0*
T0*
_output_shapes
:	�
�
Tensordot_3/MatMulMatMulTensordot_3/Reshape:output:0Tensordot_3/transpose_1:y:0*
T0*'
_output_shapes
:���������
]
Tensordot_3/Const_2Const*
_output_shapes
:*
dtype0*
valueB:
[
Tensordot_3/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concat_1ConcatV2Tensordot_3/GatherV2:output:0Tensordot_3/Const_2:output:0"Tensordot_3/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_3ReshapeTensordot_3/MatMul:product:0Tensordot_3/concat_1:output:0*
T0*'
_output_shapes
:���������
X
Mul_1Mulmul_1_xconcat:output:0*
T0*'
_output_shapes
:���������
_
addAddV2Tensordot_3:output:0	Mul_1:z:0*
T0*'
_output_shapes
:���������
f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSliceadd:z:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������	*

begin_mask*
end_mask�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977804*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977806*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977808* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977810*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977812* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977814*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977816* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977818*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977820* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977822*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977824* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977826*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977828* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977830*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977832*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
IdentityIdentity7ResNet_entropy_closure/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������i

Identity_1IdentityPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������	i

Identity_2Identitystrided_slice_2:output:0^NoOp*
T0*'
_output_shapes
:���������	�
NoOpNoOp/^ResNet_entropy_closure/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp^checked
^checked_16^dense_output/kernel/Regularizer/Square/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
2`
.ResNet_entropy_closure/StatefulPartitionedCall.ResNet_entropy_closure/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2
checkedchecked2
	checked_1	checked_12n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
�
W__inference___backward_add_1_layer_call_and_return_conditional_losses_93976666_93976683
placeholder#
gradients_add_grad_shape_inputs'
#gradients_add_grad_shape_1_inputs_1
identity

identity_1_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������g
gradients/add_grad/ShapeShapegradients_add_grad_shape_inputs*
T0*
_output_shapes
:m
gradients/add_grad/Shape_1Shape#gradients_add_grad_shape_1_inputs_1*
T0*
_output_shapes
:�
(gradients/add_grad/BroadcastGradientArgsBroadcastGradientArgs!gradients/add_grad/Shape:output:0#gradients/add_grad/Shape_1:output:0*2
_output_shapes 
:���������:����������
gradients/add_grad/SumSumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
gradients/add_grad/ReshapeReshapegradients/add_grad/Sum:output:0!gradients/add_grad/Shape:output:0*
T0*(
_output_shapes
:�����������
gradients/add_grad/Sum_1Sumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
gradients/add_grad/Reshape_1Reshape!gradients/add_grad/Sum_1:output:0#gradients/add_grad/Shape_1:output:0*
T0*(
_output_shapes
:����������l
IdentityIdentity#gradients/add_grad/Reshape:output:0*
T0*(
_output_shapes
:����������p

Identity_1Identity%gradients/add_grad/Reshape_1:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������:����������:����������*\
forward_function_nameCA__forward_add_1_layer_call_and_return_conditional_losses_93976682:. *
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������
�
�
M__inference_block_2_layer_0_layer_call_and_return_conditional_losses_93975404

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976398
input_1'
mean_shift_layer_93976249:	.
decorrelation_layer_93976252:		'
layer_input_93976255:		�#
layer_input_93976257:	�,
block_0_layer_0_93976261:
��'
block_0_layer_0_93976263:	�,
block_1_layer_0_93976268:
��'
block_1_layer_0_93976270:	�,
block_2_layer_0_93976275:
��'
block_2_layer_0_93976277:	�,
block_3_layer_0_93976282:
��'
block_3_layer_0_93976284:	�,
block_4_layer_0_93976289:
��'
block_4_layer_0_93976291:	�,
block_5_layer_0_93976296:
��'
block_5_layer_0_93976298:	�(
dense_output_93976302:	�#
dense_output_93976304:
identity��'block_0_layer_0/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_1_layer_0/StatefulPartitionedCall�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_2_layer_0/StatefulPartitionedCall�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_3_layer_0/StatefulPartitionedCall�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_4_layer_0/StatefulPartitionedCall�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_5_layer_0/StatefulPartitionedCall�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�+decorrelation_layer/StatefulPartitionedCall�$dense_output/StatefulPartitionedCall�5dense_output/kernel/Regularizer/Square/ReadVariableOp�#layer_input/StatefulPartitionedCall�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�(mean_shift_layer/StatefulPartitionedCall�
(mean_shift_layer/StatefulPartitionedCallStatefulPartitionedCallinput_1mean_shift_layer_93976249*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *W
fRRP
N__inference_mean_shift_layer_layer_call_and_return_conditional_losses_93975264�
+decorrelation_layer/StatefulPartitionedCallStatefulPartitionedCall1mean_shift_layer/StatefulPartitionedCall:output:0decorrelation_layer_93976252*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *Z
fURS
Q__inference_decorrelation_layer_layer_call_and_return_conditional_losses_93975275�
#layer_input/StatefulPartitionedCallStatefulPartitionedCall4decorrelation_layer/StatefulPartitionedCall:output:0layer_input_93976255layer_input_93976257*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *R
fMRK
I__inference_layer_input_layer_call_and_return_conditional_losses_93975301�
tf.math.softplus/SoftplusSoftplus,layer_input/StatefulPartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_0_layer_0/StatefulPartitionedCallStatefulPartitionedCall'tf.math.softplus/Softplus:activations:0block_0_layer_0_93976261block_0_layer_0_93976263*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_0_layer_0_layer_call_and_return_conditional_losses_93975330�
add/PartitionedCallPartitionedCall,layer_input/StatefulPartitionedCall:output:00block_0_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__inference_add_layer_call_and_return_conditional_losses_93975342x
tf.math.softplus_1/SoftplusSoftplusadd/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_1_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_1/Softplus:activations:0block_1_layer_0_93976268block_1_layer_0_93976270*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_1_layer_0_layer_call_and_return_conditional_losses_93975367�
add_1/PartitionedCallPartitionedCalladd/PartitionedCall:output:00block_1_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_1_layer_call_and_return_conditional_losses_93975379z
tf.math.softplus_2/SoftplusSoftplusadd_1/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_2_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_2/Softplus:activations:0block_2_layer_0_93976275block_2_layer_0_93976277*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_2_layer_0_layer_call_and_return_conditional_losses_93975404�
add_2/PartitionedCallPartitionedCalladd_1/PartitionedCall:output:00block_2_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_2_layer_call_and_return_conditional_losses_93975416z
tf.math.softplus_3/SoftplusSoftplusadd_2/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_3_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_3/Softplus:activations:0block_3_layer_0_93976282block_3_layer_0_93976284*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_3_layer_0_layer_call_and_return_conditional_losses_93975441�
add_3/PartitionedCallPartitionedCalladd_2/PartitionedCall:output:00block_3_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_3_layer_call_and_return_conditional_losses_93975453z
tf.math.softplus_4/SoftplusSoftplusadd_3/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_4_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_4/Softplus:activations:0block_4_layer_0_93976289block_4_layer_0_93976291*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_4_layer_0_layer_call_and_return_conditional_losses_93975478�
add_4/PartitionedCallPartitionedCalladd_3/PartitionedCall:output:00block_4_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_4_layer_call_and_return_conditional_losses_93975490z
tf.math.softplus_5/SoftplusSoftplusadd_4/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_5_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_5/Softplus:activations:0block_5_layer_0_93976296block_5_layer_0_93976298*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_5_layer_0_layer_call_and_return_conditional_losses_93975515�
add_5/PartitionedCallPartitionedCalladd_4/PartitionedCall:output:00block_5_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_5_layer_call_and_return_conditional_losses_93975527�
$dense_output/StatefulPartitionedCallStatefulPartitionedCalladd_5/PartitionedCall:output:0dense_output_93976302dense_output_93976304*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *S
fNRL
J__inference_dense_output_layer_call_and_return_conditional_losses_93975545�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93976255*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93976257*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93976261* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93976263*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93976268* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93976270*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93976275* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93976277*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93976282* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93976284*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93976289* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93976291*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93976296* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93976298*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_output_93976302*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: |
IdentityIdentity-dense_output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������

NoOpNoOp(^block_0_layer_0/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_1_layer_0/StatefulPartitionedCall7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_2_layer_0/StatefulPartitionedCall7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_3_layer_0/StatefulPartitionedCall7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_4_layer_0/StatefulPartitionedCall7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_5_layer_0/StatefulPartitionedCall7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp,^decorrelation_layer/StatefulPartitionedCall%^dense_output/StatefulPartitionedCall6^dense_output/kernel/Regularizer/Square/ReadVariableOp$^layer_input/StatefulPartitionedCall3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp)^mean_shift_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 2R
'block_0_layer_0/StatefulPartitionedCall'block_0_layer_0/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_1_layer_0/StatefulPartitionedCall'block_1_layer_0/StatefulPartitionedCall2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_2_layer_0/StatefulPartitionedCall'block_2_layer_0/StatefulPartitionedCall2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_3_layer_0/StatefulPartitionedCall'block_3_layer_0/StatefulPartitionedCall2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_4_layer_0/StatefulPartitionedCall'block_4_layer_0/StatefulPartitionedCall2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_5_layer_0/StatefulPartitionedCall'block_5_layer_0/StatefulPartitionedCall2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2Z
+decorrelation_layer/StatefulPartitionedCall+decorrelation_layer/StatefulPartitionedCall2L
$dense_output/StatefulPartitionedCall$dense_output/StatefulPartitionedCall2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2J
#layer_input/StatefulPartitionedCall#layer_input/StatefulPartitionedCall2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp2T
(mean_shift_layer/StatefulPartitionedCall(mean_shift_layer/StatefulPartitionedCall:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1
�
�
0__inference_sobolev_model_layer_call_fn_93978633
x
unknown:	
	unknown_0:		
	unknown_1:		�
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:
��

unknown_12:	�

unknown_13:
��

unknown_14:	�

unknown_15:	�

unknown_16:

unknown_17

unknown_18

unknown_19
identity

identity_1

identity_2��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallxunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19*!
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������	:���������	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93977124o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������	q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������	`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
22
StatefulPartitionedCallStatefulPartitionedCall:J F
'
_output_shapes
:���������	

_user_specified_namex:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
�
/__inference_dense_output_layer_call_fn_93980366

inputs
unknown:	�
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *S
fNRL
J__inference_dense_output_layer_call_and_return_conditional_losses_93975545o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
J__inference_dense_output_layer_call_and_return_conditional_losses_93975545

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�5dense_output/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: _
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp6^dense_output/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
a__inference___backward_block_1_layer_0_layer_call_and_return_conditional_losses_93976690_93976703
placeholder6
2gradients_matmul_grad_matmul_matmul_readvariableop)
%gradients_matmul_grad_matmul_1_inputs
identity

identity_1

identity_2_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������u
"gradients/BiasAdd_grad/BiasAddGradBiasAddGradgradients/grad_ys_0:output:0*
T0*
_output_shapes	
:��
gradients/MatMul_grad/MatMulMatMulgradients/grad_ys_0:output:02gradients_matmul_grad_matmul_matmul_readvariableop*
T0*(
_output_shapes
:����������*
transpose_b(�
gradients/MatMul_grad/MatMul_1MatMul%gradients_matmul_grad_matmul_1_inputsgradients/grad_ys_0:output:0*
T0* 
_output_shapes
:
��*
transpose_a(o
IdentityIdentity&gradients/MatMul_grad/MatMul:product:0*
T0*(
_output_shapes
:����������k

Identity_1Identity(gradients/MatMul_grad/MatMul_1:product:0*
T0* 
_output_shapes
:
��i

Identity_2Identity+gradients/BiasAdd_grad/BiasAddGrad:output:0*
T0*
_output_shapes	
:�"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*G
_input_shapes6
4:����������:
��:����������*f
forward_function_nameMK__forward_block_1_layer_0_layer_call_and_return_conditional_losses_93976702:. *
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������
�
T
(__inference_add_4_layer_call_fn_93980290
inputs_0
inputs_1
identity�
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_4_layer_call_and_return_conditional_losses_93975490a
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
��
�
#__inference__wrapped_model_93975250
input_1_
Qsobolev_model_resnet_entropy_closure_mean_shift_layer_sub_readvariableop_resource:	i
Wsobolev_model_resnet_entropy_closure_decorrelation_layer_matmul_readvariableop_resource:		b
Osobolev_model_resnet_entropy_closure_layer_input_matmul_readvariableop_resource:		�_
Psobolev_model_resnet_entropy_closure_layer_input_biasadd_readvariableop_resource:	�g
Ssobolev_model_resnet_entropy_closure_block_0_layer_0_matmul_readvariableop_resource:
��c
Tsobolev_model_resnet_entropy_closure_block_0_layer_0_biasadd_readvariableop_resource:	�g
Ssobolev_model_resnet_entropy_closure_block_1_layer_0_matmul_readvariableop_resource:
��c
Tsobolev_model_resnet_entropy_closure_block_1_layer_0_biasadd_readvariableop_resource:	�g
Ssobolev_model_resnet_entropy_closure_block_2_layer_0_matmul_readvariableop_resource:
��c
Tsobolev_model_resnet_entropy_closure_block_2_layer_0_biasadd_readvariableop_resource:	�g
Ssobolev_model_resnet_entropy_closure_block_3_layer_0_matmul_readvariableop_resource:
��c
Tsobolev_model_resnet_entropy_closure_block_3_layer_0_biasadd_readvariableop_resource:	�g
Ssobolev_model_resnet_entropy_closure_block_4_layer_0_matmul_readvariableop_resource:
��c
Tsobolev_model_resnet_entropy_closure_block_4_layer_0_biasadd_readvariableop_resource:	�g
Ssobolev_model_resnet_entropy_closure_block_5_layer_0_matmul_readvariableop_resource:
��c
Tsobolev_model_resnet_entropy_closure_block_5_layer_0_biasadd_readvariableop_resource:	�c
Psobolev_model_resnet_entropy_closure_dense_output_matmul_readvariableop_resource:	�_
Qsobolev_model_resnet_entropy_closure_dense_output_biasadd_readvariableop_resource:
sobolev_model_93975128
sobolev_model_tensordot_1_b
sobolev_model_mul_1_x
identity

identity_1

identity_2��Ksobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp�Jsobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp�Ksobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp�Jsobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp�Ksobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp�Jsobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp�Ksobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp�Jsobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp�Ksobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp�Jsobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp�Ksobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp�Jsobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp�Nsobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp�Hsobolev_model/ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp�Gsobolev_model/ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp�Gsobolev_model/ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp�Fsobolev_model/ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp�Hsobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp�sobolev_model/checked�sobolev_model/checked_1�
Hsobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOpReadVariableOpQsobolev_model_resnet_entropy_closure_mean_shift_layer_sub_readvariableop_resource*
_output_shapes
:	*
dtype0�
9sobolev_model/ResNet_entropy_closure/mean_shift_layer/subSubinput_1Psobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
Nsobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOpReadVariableOpWsobolev_model_resnet_entropy_closure_decorrelation_layer_matmul_readvariableop_resource*
_output_shapes

:		*
dtype0�
?sobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMulMatMul=sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub:z:0Vsobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
Fsobolev_model/ResNet_entropy_closure/layer_input/MatMul/ReadVariableOpReadVariableOpOsobolev_model_resnet_entropy_closure_layer_input_matmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
7sobolev_model/ResNet_entropy_closure/layer_input/MatMulMatMulIsobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul:product:0Nsobolev_model/ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
Gsobolev_model/ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOpReadVariableOpPsobolev_model_resnet_entropy_closure_layer_input_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
8sobolev_model/ResNet_entropy_closure/layer_input/BiasAddBiasAddAsobolev_model/ResNet_entropy_closure/layer_input/MatMul:product:0Osobolev_model/ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
>sobolev_model/ResNet_entropy_closure/tf.math.softplus/SoftplusSoftplusAsobolev_model/ResNet_entropy_closure/layer_input/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
Jsobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOpReadVariableOpSsobolev_model_resnet_entropy_closure_block_0_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
;sobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMulMatMulLsobolev_model/ResNet_entropy_closure/tf.math.softplus/Softplus:activations:0Rsobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
Ksobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOpReadVariableOpTsobolev_model_resnet_entropy_closure_block_0_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
<sobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAddBiasAddEsobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul:product:0Ssobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sobolev_model/ResNet_entropy_closure/add/addAddV2Asobolev_model/ResNet_entropy_closure/layer_input/BiasAdd:output:0Esobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
@sobolev_model/ResNet_entropy_closure/tf.math.softplus_1/SoftplusSoftplus0sobolev_model/ResNet_entropy_closure/add/add:z:0*
T0*(
_output_shapes
:�����������
Jsobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOpReadVariableOpSsobolev_model_resnet_entropy_closure_block_1_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
;sobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMulMatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_1/Softplus:activations:0Rsobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
Ksobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOpReadVariableOpTsobolev_model_resnet_entropy_closure_block_1_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
<sobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAddBiasAddEsobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul:product:0Ssobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sobolev_model/ResNet_entropy_closure/add_1/addAddV20sobolev_model/ResNet_entropy_closure/add/add:z:0Esobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
@sobolev_model/ResNet_entropy_closure/tf.math.softplus_2/SoftplusSoftplus2sobolev_model/ResNet_entropy_closure/add_1/add:z:0*
T0*(
_output_shapes
:�����������
Jsobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOpReadVariableOpSsobolev_model_resnet_entropy_closure_block_2_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
;sobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMulMatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_2/Softplus:activations:0Rsobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
Ksobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOpReadVariableOpTsobolev_model_resnet_entropy_closure_block_2_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
<sobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAddBiasAddEsobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul:product:0Ssobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sobolev_model/ResNet_entropy_closure/add_2/addAddV22sobolev_model/ResNet_entropy_closure/add_1/add:z:0Esobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
@sobolev_model/ResNet_entropy_closure/tf.math.softplus_3/SoftplusSoftplus2sobolev_model/ResNet_entropy_closure/add_2/add:z:0*
T0*(
_output_shapes
:�����������
Jsobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOpReadVariableOpSsobolev_model_resnet_entropy_closure_block_3_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
;sobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMulMatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_3/Softplus:activations:0Rsobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
Ksobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOpReadVariableOpTsobolev_model_resnet_entropy_closure_block_3_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
<sobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAddBiasAddEsobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul:product:0Ssobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sobolev_model/ResNet_entropy_closure/add_3/addAddV22sobolev_model/ResNet_entropy_closure/add_2/add:z:0Esobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
@sobolev_model/ResNet_entropy_closure/tf.math.softplus_4/SoftplusSoftplus2sobolev_model/ResNet_entropy_closure/add_3/add:z:0*
T0*(
_output_shapes
:�����������
Jsobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOpReadVariableOpSsobolev_model_resnet_entropy_closure_block_4_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
;sobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMulMatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_4/Softplus:activations:0Rsobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
Ksobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOpReadVariableOpTsobolev_model_resnet_entropy_closure_block_4_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
<sobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAddBiasAddEsobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul:product:0Ssobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sobolev_model/ResNet_entropy_closure/add_4/addAddV22sobolev_model/ResNet_entropy_closure/add_3/add:z:0Esobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
@sobolev_model/ResNet_entropy_closure/tf.math.softplus_5/SoftplusSoftplus2sobolev_model/ResNet_entropy_closure/add_4/add:z:0*
T0*(
_output_shapes
:�����������
Jsobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOpReadVariableOpSsobolev_model_resnet_entropy_closure_block_5_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
;sobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMulMatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_5/Softplus:activations:0Rsobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
Ksobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOpReadVariableOpTsobolev_model_resnet_entropy_closure_block_5_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
<sobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAddBiasAddEsobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul:product:0Ssobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sobolev_model/ResNet_entropy_closure/add_5/addAddV22sobolev_model/ResNet_entropy_closure/add_4/add:z:0Esobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
Gsobolev_model/ResNet_entropy_closure/dense_output/MatMul/ReadVariableOpReadVariableOpPsobolev_model_resnet_entropy_closure_dense_output_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
8sobolev_model/ResNet_entropy_closure/dense_output/MatMulMatMul2sobolev_model/ResNet_entropy_closure/add_5/add:z:0Osobolev_model/ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
Hsobolev_model/ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOpReadVariableOpQsobolev_model_resnet_entropy_closure_dense_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
9sobolev_model/ResNet_entropy_closure/dense_output/BiasAddBiasAddBsobolev_model/ResNet_entropy_closure/dense_output/MatMul:product:0Psobolev_model/ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
sobolev_model/ones_like/ShapeShapeBsobolev_model/ResNet_entropy_closure/dense_output/BiasAdd:output:0*
T0*
_output_shapes
:b
sobolev_model/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?�
sobolev_model/ones_likeFill&sobolev_model/ones_like/Shape:output:0&sobolev_model/ones_like/Const:output:0*
T0*'
_output_shapes
:����������
Sgradient_tape/sobolev_model/ResNet_entropy_closure/dense_output/BiasAdd/BiasAddGradBiasAddGrad sobolev_model/ones_like:output:0*
T0*
_output_shapes
:�
Mgradient_tape/sobolev_model/ResNet_entropy_closure/dense_output/MatMul/MatMulMatMul sobolev_model/ones_like:output:0Osobolev_model/ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Ogradient_tape/sobolev_model/ResNet_entropy_closure/dense_output/MatMul/MatMul_1MatMul2sobolev_model/ResNet_entropy_closure/add_5/add:z:0 sobolev_model/ones_like:output:0*
T0*
_output_shapes
:	�*
transpose_a(�
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/ShapeShape2sobolev_model/ResNet_entropy_closure/add_4/add:z:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Shape_1ShapeEsobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/BroadcastGradientArgsBroadcastGradientArgsKgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Shape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
@gradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/SumSumWgradient_tape/sobolev_model/ResNet_entropy_closure/dense_output/MatMul/MatMul:product:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/ReshapeReshapeIgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Sum:output:0Kgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Shape:output:0*
T0*(
_output_shapes
:�����������
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Sum_1SumWgradient_tape/sobolev_model/ResNet_entropy_closure/dense_output/MatMul/MatMul:product:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
Fgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Reshape_1ReshapeKgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Sum_1:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Vgradient_tape/sobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAdd/BiasAddGradBiasAddGradOgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Pgradient_tape/sobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/MatMulMatMulOgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Reshape_1:output:0Rsobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/MatMul_1MatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_5/Softplus:activations:0Ogradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_5/SigmoidSigmoid2sobolev_model/ResNet_entropy_closure/add_4/add:z:0*
T0*(
_output_shapes
:�����������
Igradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_5/mulMulZgradient_tape/sobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/MatMul:product:0Qgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_5/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
sobolev_model/AddNAddNMgradient_tape/sobolev_model/ResNet_entropy_closure/add_5/add/Reshape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_5/mul:z:0*
N*
T0*(
_output_shapes
:�����������
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/ShapeShape2sobolev_model/ResNet_entropy_closure/add_3/add:z:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Shape_1ShapeEsobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/BroadcastGradientArgsBroadcastGradientArgsKgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Shape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
@gradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/SumSumsobolev_model/AddN:sum:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/ReshapeReshapeIgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Sum:output:0Kgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Shape:output:0*
T0*(
_output_shapes
:�����������
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Sum_1Sumsobolev_model/AddN:sum:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
Fgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Reshape_1ReshapeKgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Sum_1:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Vgradient_tape/sobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAdd/BiasAddGradBiasAddGradOgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Pgradient_tape/sobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/MatMulMatMulOgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Reshape_1:output:0Rsobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/MatMul_1MatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_4/Softplus:activations:0Ogradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_4/SigmoidSigmoid2sobolev_model/ResNet_entropy_closure/add_3/add:z:0*
T0*(
_output_shapes
:�����������
Igradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_4/mulMulZgradient_tape/sobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/MatMul:product:0Qgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_4/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
sobolev_model/AddN_1AddNMgradient_tape/sobolev_model/ResNet_entropy_closure/add_4/add/Reshape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_4/mul:z:0*
N*
T0*(
_output_shapes
:�����������
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/ShapeShape2sobolev_model/ResNet_entropy_closure/add_2/add:z:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Shape_1ShapeEsobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/BroadcastGradientArgsBroadcastGradientArgsKgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Shape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
@gradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/SumSumsobolev_model/AddN_1:sum:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/ReshapeReshapeIgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Sum:output:0Kgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Shape:output:0*
T0*(
_output_shapes
:�����������
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Sum_1Sumsobolev_model/AddN_1:sum:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
Fgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Reshape_1ReshapeKgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Sum_1:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Vgradient_tape/sobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAdd/BiasAddGradBiasAddGradOgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Pgradient_tape/sobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/MatMulMatMulOgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Reshape_1:output:0Rsobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/MatMul_1MatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_3/Softplus:activations:0Ogradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_3/SigmoidSigmoid2sobolev_model/ResNet_entropy_closure/add_2/add:z:0*
T0*(
_output_shapes
:�����������
Igradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_3/mulMulZgradient_tape/sobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/MatMul:product:0Qgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_3/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
sobolev_model/AddN_2AddNMgradient_tape/sobolev_model/ResNet_entropy_closure/add_3/add/Reshape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_3/mul:z:0*
N*
T0*(
_output_shapes
:�����������
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/ShapeShape2sobolev_model/ResNet_entropy_closure/add_1/add:z:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Shape_1ShapeEsobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/BroadcastGradientArgsBroadcastGradientArgsKgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Shape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
@gradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/SumSumsobolev_model/AddN_2:sum:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/ReshapeReshapeIgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Sum:output:0Kgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Shape:output:0*
T0*(
_output_shapes
:�����������
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Sum_1Sumsobolev_model/AddN_2:sum:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
Fgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Reshape_1ReshapeKgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Sum_1:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Vgradient_tape/sobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAdd/BiasAddGradBiasAddGradOgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Pgradient_tape/sobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/MatMulMatMulOgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Reshape_1:output:0Rsobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/MatMul_1MatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_2/Softplus:activations:0Ogradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_2/SigmoidSigmoid2sobolev_model/ResNet_entropy_closure/add_1/add:z:0*
T0*(
_output_shapes
:�����������
Igradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_2/mulMulZgradient_tape/sobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/MatMul:product:0Qgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_2/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
sobolev_model/AddN_3AddNMgradient_tape/sobolev_model/ResNet_entropy_closure/add_2/add/Reshape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_2/mul:z:0*
N*
T0*(
_output_shapes
:�����������
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/ShapeShape0sobolev_model/ResNet_entropy_closure/add/add:z:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Shape_1ShapeEsobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/BroadcastGradientArgsBroadcastGradientArgsKgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Shape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
@gradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/SumSumsobolev_model/AddN_3:sum:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/ReshapeReshapeIgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Sum:output:0Kgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Shape:output:0*
T0*(
_output_shapes
:�����������
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Sum_1Sumsobolev_model/AddN_3:sum:0Wgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
Fgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Reshape_1ReshapeKgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Sum_1:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Vgradient_tape/sobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAdd/BiasAddGradBiasAddGradOgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Pgradient_tape/sobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/MatMulMatMulOgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Reshape_1:output:0Rsobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/MatMul_1MatMulNsobolev_model/ResNet_entropy_closure/tf.math.softplus_1/Softplus:activations:0Ogradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_1/SigmoidSigmoid0sobolev_model/ResNet_entropy_closure/add/add:z:0*
T0*(
_output_shapes
:�����������
Igradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_1/mulMulZgradient_tape/sobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/MatMul:product:0Qgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_1/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
sobolev_model/AddN_4AddNMgradient_tape/sobolev_model/ResNet_entropy_closure/add_1/add/Reshape:output:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus_1/mul:z:0*
N*
T0*(
_output_shapes
:�����������
@gradient_tape/sobolev_model/ResNet_entropy_closure/add/add/ShapeShapeAsobolev_model/ResNet_entropy_closure/layer_input/BiasAdd:output:0*
T0*
_output_shapes
:�
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Shape_1ShapeEsobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Pgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/BroadcastGradientArgsBroadcastGradientArgsIgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Shape:output:0Kgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
>gradient_tape/sobolev_model/ResNet_entropy_closure/add/add/SumSumsobolev_model/AddN_4:sum:0Ugradient_tape/sobolev_model/ResNet_entropy_closure/add/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
Bgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/ReshapeReshapeGgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Sum:output:0Igradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Shape:output:0*
T0*(
_output_shapes
:�����������
@gradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Sum_1Sumsobolev_model/AddN_4:sum:0Ugradient_tape/sobolev_model/ResNet_entropy_closure/add/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
Dgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Reshape_1ReshapeIgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Sum_1:output:0Kgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Vgradient_tape/sobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAdd/BiasAddGradBiasAddGradMgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Pgradient_tape/sobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/MatMulMatMulMgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Reshape_1:output:0Rsobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Rgradient_tape/sobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/MatMul_1MatMulLsobolev_model/ResNet_entropy_closure/tf.math.softplus/Softplus:activations:0Mgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
Kgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus/SigmoidSigmoidAsobolev_model/ResNet_entropy_closure/layer_input/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
Ggradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus/mulMulZgradient_tape/sobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/MatMul:product:0Ogradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
sobolev_model/AddN_5AddNKgradient_tape/sobolev_model/ResNet_entropy_closure/add/add/Reshape:output:0Kgradient_tape/sobolev_model/ResNet_entropy_closure/tf.math.softplus/mul:z:0*
N*
T0*(
_output_shapes
:�����������
Rgradient_tape/sobolev_model/ResNet_entropy_closure/layer_input/BiasAdd/BiasAddGradBiasAddGradsobolev_model/AddN_5:sum:0*
T0*
_output_shapes	
:��
Lgradient_tape/sobolev_model/ResNet_entropy_closure/layer_input/MatMul/MatMulMatMulsobolev_model/AddN_5:sum:0Nsobolev_model/ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	*
transpose_b(�
Ngradient_tape/sobolev_model/ResNet_entropy_closure/layer_input/MatMul/MatMul_1MatMulIsobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul:product:0sobolev_model/AddN_5:sum:0*
T0*
_output_shapes
:		�*
transpose_a(�
Tgradient_tape/sobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul/MatMulMatMulVgradient_tape/sobolev_model/ResNet_entropy_closure/layer_input/MatMul/MatMul:product:0Vsobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	*
transpose_b(�
Mgradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/ShapeShapeinput_1*
T0*
_output_shapes
:�
Ogradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/Shape_1ShapePsobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp:value:0*
T0*
_output_shapes
:�
]gradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/BroadcastGradientArgsBroadcastGradientArgsVgradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/Shape:output:0Xgradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/Shape_1:output:0*2
_output_shapes 
:���������:����������
Kgradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/SumSum^gradient_tape/sobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul/MatMul:product:0bgradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
Ogradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/ReshapeReshapeTgradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/Sum:output:0Vgradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/Shape:output:0*
T0*'
_output_shapes
:���������	�
sobolev_model/CastCastXgradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/Reshape:output:0*

DstT0*

SrcT0*'
_output_shapes
:���������	�
sobolev_model/checkedCheckNumericssobolev_model/Cast:y:0*
T0*'
_output_shapes
:���������	*r
messagegeinput tensor checking error at alpha = Tensor("sobolev_model/Cast:0", shape=(None, 9), dtype=float64)r
)sobolev_model/checkedandclipped/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
'sobolev_model/checkedandclipped/MinimumMinimumsobolev_model/checked:output:02sobolev_model/checkedandclipped/Minimum/y:output:0*
T0*'
_output_shapes
:���������	j
!sobolev_model/checkedandclipped/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
sobolev_model/checkedandclippedMaximum+sobolev_model/checkedandclipped/Minimum:z:0*sobolev_model/checkedandclipped/y:output:0*
T0*'
_output_shapes
:���������	r
!sobolev_model/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       t
#sobolev_model/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        t
#sobolev_model/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
sobolev_model/strided_sliceStridedSlicesobolev_model_93975128*sobolev_model/strided_slice/stack:output:0,sobolev_model/strided_slice/stack_1:output:0,sobolev_model/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:		�*

begin_mask*
end_maskf
sobolev_model/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:f
sobolev_model/Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB: p
sobolev_model/Tensordot/ShapeShape#sobolev_model/checkedandclipped:z:0*
T0*
_output_shapes
:g
%sobolev_model/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
 sobolev_model/Tensordot/GatherV2GatherV2&sobolev_model/Tensordot/Shape:output:0%sobolev_model/Tensordot/free:output:0.sobolev_model/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:i
'sobolev_model/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
"sobolev_model/Tensordot/GatherV2_1GatherV2&sobolev_model/Tensordot/Shape:output:0%sobolev_model/Tensordot/axes:output:00sobolev_model/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:g
sobolev_model/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
sobolev_model/Tensordot/ProdProd)sobolev_model/Tensordot/GatherV2:output:0&sobolev_model/Tensordot/Const:output:0*
T0*
_output_shapes
: i
sobolev_model/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
sobolev_model/Tensordot/Prod_1Prod+sobolev_model/Tensordot/GatherV2_1:output:0(sobolev_model/Tensordot/Const_1:output:0*
T0*
_output_shapes
: e
#sobolev_model/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
sobolev_model/Tensordot/concatConcatV2%sobolev_model/Tensordot/free:output:0%sobolev_model/Tensordot/axes:output:0,sobolev_model/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:�
sobolev_model/Tensordot/stackPack%sobolev_model/Tensordot/Prod:output:0'sobolev_model/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:�
!sobolev_model/Tensordot/transpose	Transpose#sobolev_model/checkedandclipped:z:0'sobolev_model/Tensordot/concat:output:0*
T0*'
_output_shapes
:���������	�
sobolev_model/Tensordot/ReshapeReshape%sobolev_model/Tensordot/transpose:y:0&sobolev_model/Tensordot/stack:output:0*
T0*0
_output_shapes
:�������������������
sobolev_model/Tensordot/MatMulMatMul(sobolev_model/Tensordot/Reshape:output:0$sobolev_model/strided_slice:output:0*
T0*(
_output_shapes
:����������j
sobolev_model/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�g
%sobolev_model/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
 sobolev_model/Tensordot/concat_1ConcatV2)sobolev_model/Tensordot/GatherV2:output:0(sobolev_model/Tensordot/Const_2:output:0.sobolev_model/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
sobolev_model/TensordotReshape(sobolev_model/Tensordot/MatMul:product:0)sobolev_model/Tensordot/concat_1:output:0*
T0*(
_output_shapes
:����������m
sobolev_model/ExpExp sobolev_model/Tensordot:output:0*
T0*(
_output_shapes
:����������h
sobolev_model/Tensordot_1/axesConst*
_output_shapes
:*
dtype0*
valueB:h
sobolev_model/Tensordot_1/freeConst*
_output_shapes
:*
dtype0*
valueB: d
sobolev_model/Tensordot_1/ShapeShapesobolev_model/Exp:y:0*
T0*
_output_shapes
:i
'sobolev_model/Tensordot_1/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
"sobolev_model/Tensordot_1/GatherV2GatherV2(sobolev_model/Tensordot_1/Shape:output:0'sobolev_model/Tensordot_1/free:output:00sobolev_model/Tensordot_1/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:k
)sobolev_model/Tensordot_1/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
$sobolev_model/Tensordot_1/GatherV2_1GatherV2(sobolev_model/Tensordot_1/Shape:output:0'sobolev_model/Tensordot_1/axes:output:02sobolev_model/Tensordot_1/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:i
sobolev_model/Tensordot_1/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
sobolev_model/Tensordot_1/ProdProd+sobolev_model/Tensordot_1/GatherV2:output:0(sobolev_model/Tensordot_1/Const:output:0*
T0*
_output_shapes
: k
!sobolev_model/Tensordot_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
 sobolev_model/Tensordot_1/Prod_1Prod-sobolev_model/Tensordot_1/GatherV2_1:output:0*sobolev_model/Tensordot_1/Const_1:output:0*
T0*
_output_shapes
: g
%sobolev_model/Tensordot_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
 sobolev_model/Tensordot_1/concatConcatV2'sobolev_model/Tensordot_1/free:output:0'sobolev_model/Tensordot_1/axes:output:0.sobolev_model/Tensordot_1/concat/axis:output:0*
N*
T0*
_output_shapes
:�
sobolev_model/Tensordot_1/stackPack'sobolev_model/Tensordot_1/Prod:output:0)sobolev_model/Tensordot_1/Prod_1:output:0*
N*
T0*
_output_shapes
:�
#sobolev_model/Tensordot_1/transpose	Transposesobolev_model/Exp:y:0)sobolev_model/Tensordot_1/concat:output:0*
T0*(
_output_shapes
:�����������
!sobolev_model/Tensordot_1/ReshapeReshape'sobolev_model/Tensordot_1/transpose:y:0(sobolev_model/Tensordot_1/stack:output:0*
T0*0
_output_shapes
:������������������{
*sobolev_model/Tensordot_1/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
%sobolev_model/Tensordot_1/transpose_1	Transposesobolev_model_tensordot_1_b3sobolev_model/Tensordot_1/transpose_1/perm:output:0*
T0*
_output_shapes
:	��
 sobolev_model/Tensordot_1/MatMulMatMul*sobolev_model/Tensordot_1/Reshape:output:0)sobolev_model/Tensordot_1/transpose_1:y:0*
T0*'
_output_shapes
:���������k
!sobolev_model/Tensordot_1/Const_2Const*
_output_shapes
:*
dtype0*
valueB:i
'sobolev_model/Tensordot_1/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
"sobolev_model/Tensordot_1/concat_1ConcatV2+sobolev_model/Tensordot_1/GatherV2:output:0*sobolev_model/Tensordot_1/Const_2:output:00sobolev_model/Tensordot_1/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
sobolev_model/Tensordot_1Reshape*sobolev_model/Tensordot_1/MatMul:product:0+sobolev_model/Tensordot_1/concat_1:output:0*
T0*'
_output_shapes
:���������n
sobolev_model/LogLog"sobolev_model/Tensordot_1:output:0*
T0*'
_output_shapes
:���������a
sobolev_model/NegNegsobolev_model/Log:y:0*
T0*'
_output_shapes
:���������[
sobolev_model/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :�
sobolev_model/concatConcatV2sobolev_model/Neg:y:0sobolev_model/Cast:y:0"sobolev_model/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������
�
sobolev_model/checked_1CheckNumericssobolev_model/concat:output:0^sobolev_model/checked*
T0*'
_output_shapes
:���������
*(
messageinput tensor checking errort
+sobolev_model/checkedandclipped_1/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
)sobolev_model/checkedandclipped_1/MinimumMinimum sobolev_model/checked_1:output:04sobolev_model/checkedandclipped_1/Minimum/y:output:0*
T0*'
_output_shapes
:���������
l
#sobolev_model/checkedandclipped_1/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
!sobolev_model/checkedandclipped_1Maximum-sobolev_model/checkedandclipped_1/Minimum:z:0,sobolev_model/checkedandclipped_1/y:output:0*
T0*'
_output_shapes
:���������
h
sobolev_model/Tensordot_2/axesConst*
_output_shapes
:*
dtype0*
valueB:h
sobolev_model/Tensordot_2/freeConst*
_output_shapes
:*
dtype0*
valueB: t
sobolev_model/Tensordot_2/ShapeShape%sobolev_model/checkedandclipped_1:z:0*
T0*
_output_shapes
:i
'sobolev_model/Tensordot_2/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
"sobolev_model/Tensordot_2/GatherV2GatherV2(sobolev_model/Tensordot_2/Shape:output:0'sobolev_model/Tensordot_2/free:output:00sobolev_model/Tensordot_2/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:k
)sobolev_model/Tensordot_2/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
$sobolev_model/Tensordot_2/GatherV2_1GatherV2(sobolev_model/Tensordot_2/Shape:output:0'sobolev_model/Tensordot_2/axes:output:02sobolev_model/Tensordot_2/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:i
sobolev_model/Tensordot_2/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
sobolev_model/Tensordot_2/ProdProd+sobolev_model/Tensordot_2/GatherV2:output:0(sobolev_model/Tensordot_2/Const:output:0*
T0*
_output_shapes
: k
!sobolev_model/Tensordot_2/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
 sobolev_model/Tensordot_2/Prod_1Prod-sobolev_model/Tensordot_2/GatherV2_1:output:0*sobolev_model/Tensordot_2/Const_1:output:0*
T0*
_output_shapes
: g
%sobolev_model/Tensordot_2/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
 sobolev_model/Tensordot_2/concatConcatV2'sobolev_model/Tensordot_2/free:output:0'sobolev_model/Tensordot_2/axes:output:0.sobolev_model/Tensordot_2/concat/axis:output:0*
N*
T0*
_output_shapes
:�
sobolev_model/Tensordot_2/stackPack'sobolev_model/Tensordot_2/Prod:output:0)sobolev_model/Tensordot_2/Prod_1:output:0*
N*
T0*
_output_shapes
:�
#sobolev_model/Tensordot_2/transpose	Transpose%sobolev_model/checkedandclipped_1:z:0)sobolev_model/Tensordot_2/concat:output:0*
T0*'
_output_shapes
:���������
�
!sobolev_model/Tensordot_2/ReshapeReshape'sobolev_model/Tensordot_2/transpose:y:0(sobolev_model/Tensordot_2/stack:output:0*
T0*0
_output_shapes
:�������������������
 sobolev_model/Tensordot_2/MatMulMatMul*sobolev_model/Tensordot_2/Reshape:output:0sobolev_model_93975128*
T0*(
_output_shapes
:����������l
!sobolev_model/Tensordot_2/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�i
'sobolev_model/Tensordot_2/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
"sobolev_model/Tensordot_2/concat_1ConcatV2+sobolev_model/Tensordot_2/GatherV2:output:0*sobolev_model/Tensordot_2/Const_2:output:00sobolev_model/Tensordot_2/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
sobolev_model/Tensordot_2Reshape*sobolev_model/Tensordot_2/MatMul:product:0+sobolev_model/Tensordot_2/concat_1:output:0*
T0*(
_output_shapes
:����������q
sobolev_model/Exp_1Exp"sobolev_model/Tensordot_2:output:0*
T0*(
_output_shapes
:�����������
sobolev_model/MulMulsobolev_model/Exp_1:y:0sobolev_model_tensordot_1_b*
T0*(
_output_shapes
:����������t
#sobolev_model/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        v
%sobolev_model/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        v
%sobolev_model/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
sobolev_model/strided_slice_1StridedSlicesobolev_model_93975128,sobolev_model/strided_slice_1/stack:output:0.sobolev_model/strided_slice_1/stack_1:output:0.sobolev_model/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
:	
�*

begin_mask*
end_maskh
sobolev_model/Tensordot_3/axesConst*
_output_shapes
:*
dtype0*
valueB:h
sobolev_model/Tensordot_3/freeConst*
_output_shapes
:*
dtype0*
valueB: d
sobolev_model/Tensordot_3/ShapeShapesobolev_model/Mul:z:0*
T0*
_output_shapes
:i
'sobolev_model/Tensordot_3/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
"sobolev_model/Tensordot_3/GatherV2GatherV2(sobolev_model/Tensordot_3/Shape:output:0'sobolev_model/Tensordot_3/free:output:00sobolev_model/Tensordot_3/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:k
)sobolev_model/Tensordot_3/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
$sobolev_model/Tensordot_3/GatherV2_1GatherV2(sobolev_model/Tensordot_3/Shape:output:0'sobolev_model/Tensordot_3/axes:output:02sobolev_model/Tensordot_3/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:i
sobolev_model/Tensordot_3/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
sobolev_model/Tensordot_3/ProdProd+sobolev_model/Tensordot_3/GatherV2:output:0(sobolev_model/Tensordot_3/Const:output:0*
T0*
_output_shapes
: k
!sobolev_model/Tensordot_3/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
 sobolev_model/Tensordot_3/Prod_1Prod-sobolev_model/Tensordot_3/GatherV2_1:output:0*sobolev_model/Tensordot_3/Const_1:output:0*
T0*
_output_shapes
: g
%sobolev_model/Tensordot_3/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
 sobolev_model/Tensordot_3/concatConcatV2'sobolev_model/Tensordot_3/free:output:0'sobolev_model/Tensordot_3/axes:output:0.sobolev_model/Tensordot_3/concat/axis:output:0*
N*
T0*
_output_shapes
:�
sobolev_model/Tensordot_3/stackPack'sobolev_model/Tensordot_3/Prod:output:0)sobolev_model/Tensordot_3/Prod_1:output:0*
N*
T0*
_output_shapes
:�
#sobolev_model/Tensordot_3/transpose	Transposesobolev_model/Mul:z:0)sobolev_model/Tensordot_3/concat:output:0*
T0*(
_output_shapes
:�����������
!sobolev_model/Tensordot_3/ReshapeReshape'sobolev_model/Tensordot_3/transpose:y:0(sobolev_model/Tensordot_3/stack:output:0*
T0*0
_output_shapes
:������������������{
*sobolev_model/Tensordot_3/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
%sobolev_model/Tensordot_3/transpose_1	Transpose&sobolev_model/strided_slice_1:output:03sobolev_model/Tensordot_3/transpose_1/perm:output:0*
T0*
_output_shapes
:	�
�
 sobolev_model/Tensordot_3/MatMulMatMul*sobolev_model/Tensordot_3/Reshape:output:0)sobolev_model/Tensordot_3/transpose_1:y:0*
T0*'
_output_shapes
:���������
k
!sobolev_model/Tensordot_3/Const_2Const*
_output_shapes
:*
dtype0*
valueB:
i
'sobolev_model/Tensordot_3/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
"sobolev_model/Tensordot_3/concat_1ConcatV2+sobolev_model/Tensordot_3/GatherV2:output:0*sobolev_model/Tensordot_3/Const_2:output:00sobolev_model/Tensordot_3/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
sobolev_model/Tensordot_3Reshape*sobolev_model/Tensordot_3/MatMul:product:0+sobolev_model/Tensordot_3/concat_1:output:0*
T0*'
_output_shapes
:���������
�
sobolev_model/Mul_1Mulsobolev_model_mul_1_xsobolev_model/concat:output:0*
T0*'
_output_shapes
:���������
�
sobolev_model/addAddV2"sobolev_model/Tensordot_3:output:0sobolev_model/Mul_1:z:0*
T0*'
_output_shapes
:���������
t
#sobolev_model/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       v
%sobolev_model/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        v
%sobolev_model/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
sobolev_model/strided_slice_2StridedSlicesobolev_model/add:z:0,sobolev_model/strided_slice_2/stack:output:0.sobolev_model/strided_slice_2/stack_1:output:0.sobolev_model/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������	*

begin_mask*
end_mask�
IdentityIdentityBsobolev_model/ResNet_entropy_closure/dense_output/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������

Identity_1IdentityXgradient_tape/sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/Reshape:output:0^NoOp*
T0*'
_output_shapes
:���������	w

Identity_2Identity&sobolev_model/strided_slice_2:output:0^NoOp*
T0*'
_output_shapes
:���������	�
NoOpNoOpL^sobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOpK^sobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOpL^sobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOpK^sobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOpL^sobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOpK^sobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOpL^sobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOpK^sobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOpL^sobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOpK^sobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOpL^sobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOpK^sobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOpO^sobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOpI^sobolev_model/ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOpH^sobolev_model/ResNet_entropy_closure/dense_output/MatMul/ReadVariableOpH^sobolev_model/ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOpG^sobolev_model/ResNet_entropy_closure/layer_input/MatMul/ReadVariableOpI^sobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp^sobolev_model/checked^sobolev_model/checked_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
2�
Ksobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOpKsobolev_model/ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp2�
Jsobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOpJsobolev_model/ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp2�
Ksobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOpKsobolev_model/ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp2�
Jsobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOpJsobolev_model/ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp2�
Ksobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOpKsobolev_model/ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp2�
Jsobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOpJsobolev_model/ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp2�
Ksobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOpKsobolev_model/ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp2�
Jsobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOpJsobolev_model/ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp2�
Ksobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOpKsobolev_model/ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp2�
Jsobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOpJsobolev_model/ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp2�
Ksobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOpKsobolev_model/ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp2�
Jsobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOpJsobolev_model/ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp2�
Nsobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOpNsobolev_model/ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp2�
Hsobolev_model/ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOpHsobolev_model/ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp2�
Gsobolev_model/ResNet_entropy_closure/dense_output/MatMul/ReadVariableOpGsobolev_model/ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp2�
Gsobolev_model/ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOpGsobolev_model/ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp2�
Fsobolev_model/ResNet_entropy_closure/layer_input/MatMul/ReadVariableOpFsobolev_model/ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp2�
Hsobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOpHsobolev_model/ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp2.
sobolev_model/checkedsobolev_model/checked22
sobolev_model/checked_1sobolev_model/checked_1:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
�
__inference_loss_fn_10_93980503U
Ablock_4_layer_0_kernel_regularizer_square_readvariableop_resource:
��
identity��8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpAblock_4_layer_0_kernel_regularizer_square_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: h
IdentityIdentity*block_4_layer_0/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp
�
�
__inference_loss_fn_13_93980536N
?block_5_layer_0_bias_regularizer_square_readvariableop_resource:	�
identity��6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp?block_5_layer_0_bias_regularizer_square_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: f
IdentityIdentity(block_5_layer_0/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 
NoOpNoOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp
�
�
2__inference_block_4_layer_0_layer_call_fn_93980262

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_4_layer_0_layer_call_and_return_conditional_losses_93975478p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
M__inference_block_1_layer_0_layer_call_and_return_conditional_losses_93980119

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93977697
x-
resnet_entropy_closure_93977229:	1
resnet_entropy_closure_93977231:		2
resnet_entropy_closure_93977233:		�.
resnet_entropy_closure_93977235:	�3
resnet_entropy_closure_93977237:
��.
resnet_entropy_closure_93977239:	�3
resnet_entropy_closure_93977241:
��.
resnet_entropy_closure_93977243:	�3
resnet_entropy_closure_93977245:
��.
resnet_entropy_closure_93977247:	�3
resnet_entropy_closure_93977249:
��.
resnet_entropy_closure_93977251:	�3
resnet_entropy_closure_93977253:
��.
resnet_entropy_closure_93977255:	�3
resnet_entropy_closure_93977257:
��.
resnet_entropy_closure_93977259:	�2
resnet_entropy_closure_93977261:	�-
resnet_entropy_closure_93977263:
unknown
tensordot_1_b
mul_1_x
identity

identity_1

identity_2��.ResNet_entropy_closure/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�checked�	checked_1�5dense_output/kernel/Regularizer/Square/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�
.ResNet_entropy_closure/StatefulPartitionedCallStatefulPartitionedCallxresnet_entropy_closure_93977229resnet_entropy_closure_93977231resnet_entropy_closure_93977233resnet_entropy_closure_93977235resnet_entropy_closure_93977237resnet_entropy_closure_93977239resnet_entropy_closure_93977241resnet_entropy_closure_93977243resnet_entropy_closure_93977245resnet_entropy_closure_93977247resnet_entropy_closure_93977249resnet_entropy_closure_93977251resnet_entropy_closure_93977253resnet_entropy_closure_93977255resnet_entropy_closure_93977257resnet_entropy_closure_93977259resnet_entropy_closure_93977261resnet_entropy_closure_93977263*
Tin
2*3
Tout+
)2'*
_collective_manager_ids
 *�
_output_shapes�
�:���������:	�:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:		�:���������	:		:���������	:���������	:	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *[
fVRT
R__forward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93977420v
ones_like/ShapeShape7ResNet_entropy_closure/StatefulPartitionedCall:output:0*
T0*
_output_shapes
:T
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?w
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:����������
PartitionedCallPartitionedCallones_like:output:07ResNet_entropy_closure/StatefulPartitionedCall:output:17ResNet_entropy_closure/StatefulPartitionedCall:output:27ResNet_entropy_closure/StatefulPartitionedCall:output:37ResNet_entropy_closure/StatefulPartitionedCall:output:47ResNet_entropy_closure/StatefulPartitionedCall:output:57ResNet_entropy_closure/StatefulPartitionedCall:output:67ResNet_entropy_closure/StatefulPartitionedCall:output:77ResNet_entropy_closure/StatefulPartitionedCall:output:87ResNet_entropy_closure/StatefulPartitionedCall:output:98ResNet_entropy_closure/StatefulPartitionedCall:output:108ResNet_entropy_closure/StatefulPartitionedCall:output:118ResNet_entropy_closure/StatefulPartitionedCall:output:128ResNet_entropy_closure/StatefulPartitionedCall:output:138ResNet_entropy_closure/StatefulPartitionedCall:output:148ResNet_entropy_closure/StatefulPartitionedCall:output:158ResNet_entropy_closure/StatefulPartitionedCall:output:168ResNet_entropy_closure/StatefulPartitionedCall:output:178ResNet_entropy_closure/StatefulPartitionedCall:output:188ResNet_entropy_closure/StatefulPartitionedCall:output:198ResNet_entropy_closure/StatefulPartitionedCall:output:208ResNet_entropy_closure/StatefulPartitionedCall:output:218ResNet_entropy_closure/StatefulPartitionedCall:output:228ResNet_entropy_closure/StatefulPartitionedCall:output:238ResNet_entropy_closure/StatefulPartitionedCall:output:248ResNet_entropy_closure/StatefulPartitionedCall:output:258ResNet_entropy_closure/StatefulPartitionedCall:output:268ResNet_entropy_closure/StatefulPartitionedCall:output:278ResNet_entropy_closure/StatefulPartitionedCall:output:288ResNet_entropy_closure/StatefulPartitionedCall:output:298ResNet_entropy_closure/StatefulPartitionedCall:output:308ResNet_entropy_closure/StatefulPartitionedCall:output:318ResNet_entropy_closure/StatefulPartitionedCall:output:328ResNet_entropy_closure/StatefulPartitionedCall:output:338ResNet_entropy_closure/StatefulPartitionedCall:output:348ResNet_entropy_closure/StatefulPartitionedCall:output:358ResNet_entropy_closure/StatefulPartitionedCall:output:368ResNet_entropy_closure/StatefulPartitionedCall:output:378ResNet_entropy_closure/StatefulPartitionedCall:output:38*2
Tin+
)2'*
Tout
2*
_collective_manager_ids
 *�
_output_shapes�
�:���������	:	:		:		�:�:
��:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�:* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *q
flRj
h__inference___backward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93977269_93977421g
CastCastPartitionedCall:output:0*

DstT0*

SrcT0*'
_output_shapes
:���������	�
checkedCheckNumericsCast:y:0*
T0*'
_output_shapes
:���������	*d
messageYWinput tensor checking error at alpha = Tensor("Cast:0", shape=(None, 9), dtype=float64)d
checkedandclipped/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped/MinimumMinimumchecked:output:0$checkedandclipped/Minimum/y:output:0*
T0*'
_output_shapes
:���������	\
checkedandclipped/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclippedMaximumcheckedandclipped/Minimum:z:0checkedandclipped/y:output:0*
T0*'
_output_shapes
:���������	d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSliceunknownstrided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:		�*

begin_mask*
end_maskX
Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:X
Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB: T
Tensordot/ShapeShapecheckedandclipped:z:0*
T0*
_output_shapes
:Y
Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2GatherV2Tensordot/Shape:output:0Tensordot/free:output:0 Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2_1GatherV2Tensordot/Shape:output:0Tensordot/axes:output:0"Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:Y
Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: n
Tensordot/ProdProdTensordot/GatherV2:output:0Tensordot/Const:output:0*
T0*
_output_shapes
: [
Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: t
Tensordot/Prod_1ProdTensordot/GatherV2_1:output:0Tensordot/Const_1:output:0*
T0*
_output_shapes
: W
Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concatConcatV2Tensordot/free:output:0Tensordot/axes:output:0Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:y
Tensordot/stackPackTensordot/Prod:output:0Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot/transpose	Transposecheckedandclipped:z:0Tensordot/concat:output:0*
T0*'
_output_shapes
:���������	�
Tensordot/ReshapeReshapeTensordot/transpose:y:0Tensordot/stack:output:0*
T0*0
_output_shapes
:�������������������
Tensordot/MatMulMatMulTensordot/Reshape:output:0strided_slice:output:0*
T0*(
_output_shapes
:����������\
Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�Y
Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concat_1ConcatV2Tensordot/GatherV2:output:0Tensordot/Const_2:output:0 Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
	TensordotReshapeTensordot/MatMul:product:0Tensordot/concat_1:output:0*
T0*(
_output_shapes
:����������Q
ExpExpTensordot:output:0*
T0*(
_output_shapes
:����������Z
Tensordot_1/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_1/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_1/ShapeShapeExp:y:0*
T0*
_output_shapes
:[
Tensordot_1/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2GatherV2Tensordot_1/Shape:output:0Tensordot_1/free:output:0"Tensordot_1/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_1/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2_1GatherV2Tensordot_1/Shape:output:0Tensordot_1/axes:output:0$Tensordot_1/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_1/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_1/ProdProdTensordot_1/GatherV2:output:0Tensordot_1/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_1/Prod_1ProdTensordot_1/GatherV2_1:output:0Tensordot_1/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concatConcatV2Tensordot_1/free:output:0Tensordot_1/axes:output:0 Tensordot_1/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_1/stackPackTensordot_1/Prod:output:0Tensordot_1/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_1/transpose	TransposeExp:y:0Tensordot_1/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_1/ReshapeReshapeTensordot_1/transpose:y:0Tensordot_1/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_1/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_1/transpose_1	Transposetensordot_1_b%Tensordot_1/transpose_1/perm:output:0*
T0*
_output_shapes
:	��
Tensordot_1/MatMulMatMulTensordot_1/Reshape:output:0Tensordot_1/transpose_1:y:0*
T0*'
_output_shapes
:���������]
Tensordot_1/Const_2Const*
_output_shapes
:*
dtype0*
valueB:[
Tensordot_1/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concat_1ConcatV2Tensordot_1/GatherV2:output:0Tensordot_1/Const_2:output:0"Tensordot_1/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_1ReshapeTensordot_1/MatMul:product:0Tensordot_1/concat_1:output:0*
T0*'
_output_shapes
:���������R
LogLogTensordot_1:output:0*
T0*'
_output_shapes
:���������E
NegNegLog:y:0*
T0*'
_output_shapes
:���������M
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :v
concatConcatV2Neg:y:0Cast:y:0concat/axis:output:0*
N*
T0*'
_output_shapes
:���������
�
	checked_1CheckNumericsconcat:output:0^checked*
T0*'
_output_shapes
:���������
*(
messageinput tensor checking errorf
checkedandclipped_1/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped_1/MinimumMinimumchecked_1:output:0&checkedandclipped_1/Minimum/y:output:0*
T0*'
_output_shapes
:���������
^
checkedandclipped_1/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclipped_1Maximumcheckedandclipped_1/Minimum:z:0checkedandclipped_1/y:output:0*
T0*'
_output_shapes
:���������
Z
Tensordot_2/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_2/freeConst*
_output_shapes
:*
dtype0*
valueB: X
Tensordot_2/ShapeShapecheckedandclipped_1:z:0*
T0*
_output_shapes
:[
Tensordot_2/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2GatherV2Tensordot_2/Shape:output:0Tensordot_2/free:output:0"Tensordot_2/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_2/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2_1GatherV2Tensordot_2/Shape:output:0Tensordot_2/axes:output:0$Tensordot_2/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_2/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_2/ProdProdTensordot_2/GatherV2:output:0Tensordot_2/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_2/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_2/Prod_1ProdTensordot_2/GatherV2_1:output:0Tensordot_2/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_2/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concatConcatV2Tensordot_2/free:output:0Tensordot_2/axes:output:0 Tensordot_2/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_2/stackPackTensordot_2/Prod:output:0Tensordot_2/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2/transpose	Transposecheckedandclipped_1:z:0Tensordot_2/concat:output:0*
T0*'
_output_shapes
:���������
�
Tensordot_2/ReshapeReshapeTensordot_2/transpose:y:0Tensordot_2/stack:output:0*
T0*0
_output_shapes
:������������������v
Tensordot_2/MatMulMatMulTensordot_2/Reshape:output:0unknown*
T0*(
_output_shapes
:����������^
Tensordot_2/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�[
Tensordot_2/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concat_1ConcatV2Tensordot_2/GatherV2:output:0Tensordot_2/Const_2:output:0"Tensordot_2/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2ReshapeTensordot_2/MatMul:product:0Tensordot_2/concat_1:output:0*
T0*(
_output_shapes
:����������U
Exp_1ExpTensordot_2:output:0*
T0*(
_output_shapes
:����������W
MulMul	Exp_1:y:0tensordot_1_b*
T0*(
_output_shapes
:����������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSliceunknownstrided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
:	
�*

begin_mask*
end_maskZ
Tensordot_3/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_3/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_3/ShapeShapeMul:z:0*
T0*
_output_shapes
:[
Tensordot_3/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2GatherV2Tensordot_3/Shape:output:0Tensordot_3/free:output:0"Tensordot_3/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_3/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2_1GatherV2Tensordot_3/Shape:output:0Tensordot_3/axes:output:0$Tensordot_3/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_3/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_3/ProdProdTensordot_3/GatherV2:output:0Tensordot_3/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_3/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_3/Prod_1ProdTensordot_3/GatherV2_1:output:0Tensordot_3/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_3/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concatConcatV2Tensordot_3/free:output:0Tensordot_3/axes:output:0 Tensordot_3/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_3/stackPackTensordot_3/Prod:output:0Tensordot_3/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_3/transpose	TransposeMul:z:0Tensordot_3/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_3/ReshapeReshapeTensordot_3/transpose:y:0Tensordot_3/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_3/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_3/transpose_1	Transposestrided_slice_1:output:0%Tensordot_3/transpose_1/perm:output:0*
T0*
_output_shapes
:	�
�
Tensordot_3/MatMulMatMulTensordot_3/Reshape:output:0Tensordot_3/transpose_1:y:0*
T0*'
_output_shapes
:���������
]
Tensordot_3/Const_2Const*
_output_shapes
:*
dtype0*
valueB:
[
Tensordot_3/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concat_1ConcatV2Tensordot_3/GatherV2:output:0Tensordot_3/Const_2:output:0"Tensordot_3/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_3ReshapeTensordot_3/MatMul:product:0Tensordot_3/concat_1:output:0*
T0*'
_output_shapes
:���������
X
Mul_1Mulmul_1_xconcat:output:0*
T0*'
_output_shapes
:���������
_
addAddV2Tensordot_3:output:0	Mul_1:z:0*
T0*'
_output_shapes
:���������
f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSliceadd:z:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������	*

begin_mask*
end_mask�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977233*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977235*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977237* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977239*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977241* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977243*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977245* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977247*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977249* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977251*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977253* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977255*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977257* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977259*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93977261*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
IdentityIdentity7ResNet_entropy_closure/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������i

Identity_1IdentityPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������	i

Identity_2Identitystrided_slice_2:output:0^NoOp*
T0*'
_output_shapes
:���������	�
NoOpNoOp/^ResNet_entropy_closure/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp^checked
^checked_16^dense_output/kernel/Regularizer/Square/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
2`
.ResNet_entropy_closure/StatefulPartitionedCall.ResNet_entropy_closure/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2
checkedchecked2
	checked_1	checked_12n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp:J F
'
_output_shapes
:���������	

_user_specified_namex:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�

�
e__inference___backward_decorrelation_layer_layer_call_and_return_conditional_losses_93976785_93976796
placeholder6
2gradients_matmul_grad_matmul_matmul_readvariableop)
%gradients_matmul_grad_matmul_1_inputs
identity

identity_1^
gradients/grad_ys_0Identityplaceholder*
T0*'
_output_shapes
:���������	�
gradients/MatMul_grad/MatMulMatMulgradients/grad_ys_0:output:02gradients_matmul_grad_matmul_matmul_readvariableop*
T0*'
_output_shapes
:���������	*
transpose_b(�
gradients/MatMul_grad/MatMul_1MatMul%gradients_matmul_grad_matmul_1_inputsgradients/grad_ys_0:output:0*
T0*
_output_shapes

:		*
transpose_a(n
IdentityIdentity&gradients/MatMul_grad/MatMul:product:0*
T0*'
_output_shapes
:���������	i

Identity_1Identity(gradients/MatMul_grad/MatMul_1:product:0*
T0*
_output_shapes

:		"
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:���������	:		:���������	*j
forward_function_nameQO__forward_decorrelation_layer_layer_call_and_return_conditional_losses_93976795:- )
'
_output_shapes
:���������	:$ 

_output_shapes

:		:-)
'
_output_shapes
:���������	
�
�
W__inference___backward_add_4_layer_call_and_return_conditional_losses_93976519_93976536
placeholder#
gradients_add_grad_shape_inputs'
#gradients_add_grad_shape_1_inputs_1
identity

identity_1_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������g
gradients/add_grad/ShapeShapegradients_add_grad_shape_inputs*
T0*
_output_shapes
:m
gradients/add_grad/Shape_1Shape#gradients_add_grad_shape_1_inputs_1*
T0*
_output_shapes
:�
(gradients/add_grad/BroadcastGradientArgsBroadcastGradientArgs!gradients/add_grad/Shape:output:0#gradients/add_grad/Shape_1:output:0*2
_output_shapes 
:���������:����������
gradients/add_grad/SumSumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
gradients/add_grad/ReshapeReshapegradients/add_grad/Sum:output:0!gradients/add_grad/Shape:output:0*
T0*(
_output_shapes
:�����������
gradients/add_grad/Sum_1Sumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
gradients/add_grad/Reshape_1Reshape!gradients/add_grad/Sum_1:output:0#gradients/add_grad/Shape_1:output:0*
T0*(
_output_shapes
:����������l
IdentityIdentity#gradients/add_grad/Reshape:output:0*
T0*(
_output_shapes
:����������p

Identity_1Identity%gradients/add_grad/Reshape_1:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������:����������:����������*\
forward_function_nameCA__forward_add_4_layer_call_and_return_conditional_losses_93976535:. *
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������
�
�
Q__inference_decorrelation_layer_layer_call_and_return_conditional_losses_93979978

inputs0
matmul_readvariableop_resource:		
identity��MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:		*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	_
IdentityIdentityMatMul:product:0^NoOp*
T0*'
_output_shapes
:���������	^
NoOpNoOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*(
_input_shapes
:���������	: 2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
A__forward_add_1_layer_call_and_return_conditional_losses_93976682
inputs_0

inputs_1_0
identity

inputs
inputs_1U
addAddV2inputs_0
inputs_1_0*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"
inputsinputs_0"
inputs_1
inputs_1_0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������*s
backward_function_nameYW__inference___backward_add_1_layer_call_and_return_conditional_losses_93976666_93976683:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93978433
input_1-
resnet_entropy_closure_93978118:	1
resnet_entropy_closure_93978120:		2
resnet_entropy_closure_93978122:		�.
resnet_entropy_closure_93978124:	�3
resnet_entropy_closure_93978126:
��.
resnet_entropy_closure_93978128:	�3
resnet_entropy_closure_93978130:
��.
resnet_entropy_closure_93978132:	�3
resnet_entropy_closure_93978134:
��.
resnet_entropy_closure_93978136:	�3
resnet_entropy_closure_93978138:
��.
resnet_entropy_closure_93978140:	�3
resnet_entropy_closure_93978142:
��.
resnet_entropy_closure_93978144:	�3
resnet_entropy_closure_93978146:
��.
resnet_entropy_closure_93978148:	�2
resnet_entropy_closure_93978150:	�-
resnet_entropy_closure_93978152:
unknown
tensordot_1_b
mul_1_x
identity

identity_1

identity_2��.ResNet_entropy_closure/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�checked�	checked_1�5dense_output/kernel/Regularizer/Square/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�
.ResNet_entropy_closure/StatefulPartitionedCallStatefulPartitionedCallinput_1resnet_entropy_closure_93978118resnet_entropy_closure_93978120resnet_entropy_closure_93978122resnet_entropy_closure_93978124resnet_entropy_closure_93978126resnet_entropy_closure_93978128resnet_entropy_closure_93978130resnet_entropy_closure_93978132resnet_entropy_closure_93978134resnet_entropy_closure_93978136resnet_entropy_closure_93978138resnet_entropy_closure_93978140resnet_entropy_closure_93978142resnet_entropy_closure_93978144resnet_entropy_closure_93978146resnet_entropy_closure_93978148resnet_entropy_closure_93978150resnet_entropy_closure_93978152*
Tin
2*3
Tout+
)2'*
_collective_manager_ids
 *�
_output_shapes�
�:���������:	�:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:		�:���������	:		:���������	:���������	:	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *[
fVRT
R__forward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93977420v
ones_like/ShapeShape7ResNet_entropy_closure/StatefulPartitionedCall:output:0*
T0*
_output_shapes
:T
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?w
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:����������
PartitionedCallPartitionedCallones_like:output:07ResNet_entropy_closure/StatefulPartitionedCall:output:17ResNet_entropy_closure/StatefulPartitionedCall:output:27ResNet_entropy_closure/StatefulPartitionedCall:output:37ResNet_entropy_closure/StatefulPartitionedCall:output:47ResNet_entropy_closure/StatefulPartitionedCall:output:57ResNet_entropy_closure/StatefulPartitionedCall:output:67ResNet_entropy_closure/StatefulPartitionedCall:output:77ResNet_entropy_closure/StatefulPartitionedCall:output:87ResNet_entropy_closure/StatefulPartitionedCall:output:98ResNet_entropy_closure/StatefulPartitionedCall:output:108ResNet_entropy_closure/StatefulPartitionedCall:output:118ResNet_entropy_closure/StatefulPartitionedCall:output:128ResNet_entropy_closure/StatefulPartitionedCall:output:138ResNet_entropy_closure/StatefulPartitionedCall:output:148ResNet_entropy_closure/StatefulPartitionedCall:output:158ResNet_entropy_closure/StatefulPartitionedCall:output:168ResNet_entropy_closure/StatefulPartitionedCall:output:178ResNet_entropy_closure/StatefulPartitionedCall:output:188ResNet_entropy_closure/StatefulPartitionedCall:output:198ResNet_entropy_closure/StatefulPartitionedCall:output:208ResNet_entropy_closure/StatefulPartitionedCall:output:218ResNet_entropy_closure/StatefulPartitionedCall:output:228ResNet_entropy_closure/StatefulPartitionedCall:output:238ResNet_entropy_closure/StatefulPartitionedCall:output:248ResNet_entropy_closure/StatefulPartitionedCall:output:258ResNet_entropy_closure/StatefulPartitionedCall:output:268ResNet_entropy_closure/StatefulPartitionedCall:output:278ResNet_entropy_closure/StatefulPartitionedCall:output:288ResNet_entropy_closure/StatefulPartitionedCall:output:298ResNet_entropy_closure/StatefulPartitionedCall:output:308ResNet_entropy_closure/StatefulPartitionedCall:output:318ResNet_entropy_closure/StatefulPartitionedCall:output:328ResNet_entropy_closure/StatefulPartitionedCall:output:338ResNet_entropy_closure/StatefulPartitionedCall:output:348ResNet_entropy_closure/StatefulPartitionedCall:output:358ResNet_entropy_closure/StatefulPartitionedCall:output:368ResNet_entropy_closure/StatefulPartitionedCall:output:378ResNet_entropy_closure/StatefulPartitionedCall:output:38*2
Tin+
)2'*
Tout
2*
_collective_manager_ids
 *�
_output_shapes�
�:���������	:	:		:		�:�:
��:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�:* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *q
flRj
h__inference___backward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93977269_93977421g
CastCastPartitionedCall:output:0*

DstT0*

SrcT0*'
_output_shapes
:���������	�
checkedCheckNumericsCast:y:0*
T0*'
_output_shapes
:���������	*d
messageYWinput tensor checking error at alpha = Tensor("Cast:0", shape=(None, 9), dtype=float64)d
checkedandclipped/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped/MinimumMinimumchecked:output:0$checkedandclipped/Minimum/y:output:0*
T0*'
_output_shapes
:���������	\
checkedandclipped/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclippedMaximumcheckedandclipped/Minimum:z:0checkedandclipped/y:output:0*
T0*'
_output_shapes
:���������	d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSliceunknownstrided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:		�*

begin_mask*
end_maskX
Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:X
Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB: T
Tensordot/ShapeShapecheckedandclipped:z:0*
T0*
_output_shapes
:Y
Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2GatherV2Tensordot/Shape:output:0Tensordot/free:output:0 Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2_1GatherV2Tensordot/Shape:output:0Tensordot/axes:output:0"Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:Y
Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: n
Tensordot/ProdProdTensordot/GatherV2:output:0Tensordot/Const:output:0*
T0*
_output_shapes
: [
Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: t
Tensordot/Prod_1ProdTensordot/GatherV2_1:output:0Tensordot/Const_1:output:0*
T0*
_output_shapes
: W
Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concatConcatV2Tensordot/free:output:0Tensordot/axes:output:0Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:y
Tensordot/stackPackTensordot/Prod:output:0Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot/transpose	Transposecheckedandclipped:z:0Tensordot/concat:output:0*
T0*'
_output_shapes
:���������	�
Tensordot/ReshapeReshapeTensordot/transpose:y:0Tensordot/stack:output:0*
T0*0
_output_shapes
:�������������������
Tensordot/MatMulMatMulTensordot/Reshape:output:0strided_slice:output:0*
T0*(
_output_shapes
:����������\
Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�Y
Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concat_1ConcatV2Tensordot/GatherV2:output:0Tensordot/Const_2:output:0 Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
	TensordotReshapeTensordot/MatMul:product:0Tensordot/concat_1:output:0*
T0*(
_output_shapes
:����������Q
ExpExpTensordot:output:0*
T0*(
_output_shapes
:����������Z
Tensordot_1/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_1/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_1/ShapeShapeExp:y:0*
T0*
_output_shapes
:[
Tensordot_1/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2GatherV2Tensordot_1/Shape:output:0Tensordot_1/free:output:0"Tensordot_1/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_1/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2_1GatherV2Tensordot_1/Shape:output:0Tensordot_1/axes:output:0$Tensordot_1/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_1/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_1/ProdProdTensordot_1/GatherV2:output:0Tensordot_1/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_1/Prod_1ProdTensordot_1/GatherV2_1:output:0Tensordot_1/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concatConcatV2Tensordot_1/free:output:0Tensordot_1/axes:output:0 Tensordot_1/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_1/stackPackTensordot_1/Prod:output:0Tensordot_1/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_1/transpose	TransposeExp:y:0Tensordot_1/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_1/ReshapeReshapeTensordot_1/transpose:y:0Tensordot_1/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_1/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_1/transpose_1	Transposetensordot_1_b%Tensordot_1/transpose_1/perm:output:0*
T0*
_output_shapes
:	��
Tensordot_1/MatMulMatMulTensordot_1/Reshape:output:0Tensordot_1/transpose_1:y:0*
T0*'
_output_shapes
:���������]
Tensordot_1/Const_2Const*
_output_shapes
:*
dtype0*
valueB:[
Tensordot_1/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concat_1ConcatV2Tensordot_1/GatherV2:output:0Tensordot_1/Const_2:output:0"Tensordot_1/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_1ReshapeTensordot_1/MatMul:product:0Tensordot_1/concat_1:output:0*
T0*'
_output_shapes
:���������R
LogLogTensordot_1:output:0*
T0*'
_output_shapes
:���������E
NegNegLog:y:0*
T0*'
_output_shapes
:���������M
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :v
concatConcatV2Neg:y:0Cast:y:0concat/axis:output:0*
N*
T0*'
_output_shapes
:���������
�
	checked_1CheckNumericsconcat:output:0^checked*
T0*'
_output_shapes
:���������
*(
messageinput tensor checking errorf
checkedandclipped_1/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped_1/MinimumMinimumchecked_1:output:0&checkedandclipped_1/Minimum/y:output:0*
T0*'
_output_shapes
:���������
^
checkedandclipped_1/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclipped_1Maximumcheckedandclipped_1/Minimum:z:0checkedandclipped_1/y:output:0*
T0*'
_output_shapes
:���������
Z
Tensordot_2/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_2/freeConst*
_output_shapes
:*
dtype0*
valueB: X
Tensordot_2/ShapeShapecheckedandclipped_1:z:0*
T0*
_output_shapes
:[
Tensordot_2/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2GatherV2Tensordot_2/Shape:output:0Tensordot_2/free:output:0"Tensordot_2/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_2/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2_1GatherV2Tensordot_2/Shape:output:0Tensordot_2/axes:output:0$Tensordot_2/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_2/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_2/ProdProdTensordot_2/GatherV2:output:0Tensordot_2/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_2/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_2/Prod_1ProdTensordot_2/GatherV2_1:output:0Tensordot_2/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_2/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concatConcatV2Tensordot_2/free:output:0Tensordot_2/axes:output:0 Tensordot_2/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_2/stackPackTensordot_2/Prod:output:0Tensordot_2/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2/transpose	Transposecheckedandclipped_1:z:0Tensordot_2/concat:output:0*
T0*'
_output_shapes
:���������
�
Tensordot_2/ReshapeReshapeTensordot_2/transpose:y:0Tensordot_2/stack:output:0*
T0*0
_output_shapes
:������������������v
Tensordot_2/MatMulMatMulTensordot_2/Reshape:output:0unknown*
T0*(
_output_shapes
:����������^
Tensordot_2/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�[
Tensordot_2/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concat_1ConcatV2Tensordot_2/GatherV2:output:0Tensordot_2/Const_2:output:0"Tensordot_2/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2ReshapeTensordot_2/MatMul:product:0Tensordot_2/concat_1:output:0*
T0*(
_output_shapes
:����������U
Exp_1ExpTensordot_2:output:0*
T0*(
_output_shapes
:����������W
MulMul	Exp_1:y:0tensordot_1_b*
T0*(
_output_shapes
:����������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSliceunknownstrided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
:	
�*

begin_mask*
end_maskZ
Tensordot_3/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_3/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_3/ShapeShapeMul:z:0*
T0*
_output_shapes
:[
Tensordot_3/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2GatherV2Tensordot_3/Shape:output:0Tensordot_3/free:output:0"Tensordot_3/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_3/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2_1GatherV2Tensordot_3/Shape:output:0Tensordot_3/axes:output:0$Tensordot_3/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_3/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_3/ProdProdTensordot_3/GatherV2:output:0Tensordot_3/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_3/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_3/Prod_1ProdTensordot_3/GatherV2_1:output:0Tensordot_3/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_3/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concatConcatV2Tensordot_3/free:output:0Tensordot_3/axes:output:0 Tensordot_3/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_3/stackPackTensordot_3/Prod:output:0Tensordot_3/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_3/transpose	TransposeMul:z:0Tensordot_3/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_3/ReshapeReshapeTensordot_3/transpose:y:0Tensordot_3/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_3/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_3/transpose_1	Transposestrided_slice_1:output:0%Tensordot_3/transpose_1/perm:output:0*
T0*
_output_shapes
:	�
�
Tensordot_3/MatMulMatMulTensordot_3/Reshape:output:0Tensordot_3/transpose_1:y:0*
T0*'
_output_shapes
:���������
]
Tensordot_3/Const_2Const*
_output_shapes
:*
dtype0*
valueB:
[
Tensordot_3/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concat_1ConcatV2Tensordot_3/GatherV2:output:0Tensordot_3/Const_2:output:0"Tensordot_3/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_3ReshapeTensordot_3/MatMul:product:0Tensordot_3/concat_1:output:0*
T0*'
_output_shapes
:���������
X
Mul_1Mulmul_1_xconcat:output:0*
T0*'
_output_shapes
:���������
_
addAddV2Tensordot_3:output:0	Mul_1:z:0*
T0*'
_output_shapes
:���������
f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSliceadd:z:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������	*

begin_mask*
end_mask�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978122*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978124*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978126* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978128*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978130* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978132*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978134* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978136*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978138* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978140*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978142* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978144*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978146* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978148*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93978150*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
IdentityIdentity7ResNet_entropy_closure/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������i

Identity_1IdentityPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������	i

Identity_2Identitystrided_slice_2:output:0^NoOp*
T0*'
_output_shapes
:���������	�
NoOpNoOp/^ResNet_entropy_closure/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp^checked
^checked_16^dense_output/kernel/Regularizer/Square/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
2`
.ResNet_entropy_closure/StatefulPartitionedCall.ResNet_entropy_closure/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2
checkedchecked2
	checked_1	checked_12n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
�
A__forward_add_2_layer_call_and_return_conditional_losses_93976633
inputs_0

inputs_1_0
identity

inputs
inputs_1U
addAddV2inputs_0
inputs_1_0*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"
inputsinputs_0"
inputs_1
inputs_1_0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������*s
backward_function_nameYW__inference___backward_add_2_layer_call_and_return_conditional_losses_93976617_93976634:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
m
A__inference_add_layer_call_and_return_conditional_losses_93980076
inputs_0
inputs_1
identityS
addAddV2inputs_0inputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
�
�
O__forward_decorrelation_layer_layer_call_and_return_conditional_losses_93976795
inputs_00
matmul_readvariableop_resource:		
identity
matmul_readvariableop

inputs��MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:		*
dtype0k
MatMulMatMulinputs_0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	_
IdentityIdentityMatMul:product:0^NoOp*
T0*'
_output_shapes
:���������	^
NoOpNoOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"6
matmul_readvariableopMatMul/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime*(
_input_shapes
:���������	: *�
backward_function_namege__inference___backward_decorrelation_layer_layer_call_and_return_conditional_losses_93976785_939767962.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
a__inference___backward_block_4_layer_0_layer_call_and_return_conditional_losses_93976543_93976556
placeholder6
2gradients_matmul_grad_matmul_matmul_readvariableop)
%gradients_matmul_grad_matmul_1_inputs
identity

identity_1

identity_2_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������u
"gradients/BiasAdd_grad/BiasAddGradBiasAddGradgradients/grad_ys_0:output:0*
T0*
_output_shapes	
:��
gradients/MatMul_grad/MatMulMatMulgradients/grad_ys_0:output:02gradients_matmul_grad_matmul_matmul_readvariableop*
T0*(
_output_shapes
:����������*
transpose_b(�
gradients/MatMul_grad/MatMul_1MatMul%gradients_matmul_grad_matmul_1_inputsgradients/grad_ys_0:output:0*
T0* 
_output_shapes
:
��*
transpose_a(o
IdentityIdentity&gradients/MatMul_grad/MatMul:product:0*
T0*(
_output_shapes
:����������k

Identity_1Identity(gradients/MatMul_grad/MatMul_1:product:0*
T0* 
_output_shapes
:
��i

Identity_2Identity+gradients/BiasAdd_grad/BiasAddGrad:output:0*
T0*
_output_shapes	
:�"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*G
_input_shapes6
4:����������:
��:����������*f
forward_function_nameMK__forward_block_4_layer_0_layer_call_and_return_conditional_losses_93976555:. *
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������
��
�
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93979071
xQ
Cresnet_entropy_closure_mean_shift_layer_sub_readvariableop_resource:	[
Iresnet_entropy_closure_decorrelation_layer_matmul_readvariableop_resource:		T
Aresnet_entropy_closure_layer_input_matmul_readvariableop_resource:		�Q
Bresnet_entropy_closure_layer_input_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_0_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_0_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_1_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_1_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_2_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_2_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_3_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_3_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_4_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_4_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_5_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_5_layer_0_biasadd_readvariableop_resource:	�U
Bresnet_entropy_closure_dense_output_matmul_readvariableop_resource:	�Q
Cresnet_entropy_closure_dense_output_biasadd_readvariableop_resource:
unknown
tensordot_1_b
mul_1_x
identity

identity_1

identity_2��=ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp�@ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp�:ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp�9ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp�9ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp�8ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp�:ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�checked�	checked_1�5dense_output/kernel/Regularizer/Square/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�
:ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOpReadVariableOpCresnet_entropy_closure_mean_shift_layer_sub_readvariableop_resource*
_output_shapes
:	*
dtype0�
+ResNet_entropy_closure/mean_shift_layer/subSubxBResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
@ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOpReadVariableOpIresnet_entropy_closure_decorrelation_layer_matmul_readvariableop_resource*
_output_shapes

:		*
dtype0�
1ResNet_entropy_closure/decorrelation_layer/MatMulMatMul/ResNet_entropy_closure/mean_shift_layer/sub:z:0HResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
8ResNet_entropy_closure/layer_input/MatMul/ReadVariableOpReadVariableOpAresnet_entropy_closure_layer_input_matmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
)ResNet_entropy_closure/layer_input/MatMulMatMul;ResNet_entropy_closure/decorrelation_layer/MatMul:product:0@ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
9ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOpReadVariableOpBresnet_entropy_closure_layer_input_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
*ResNet_entropy_closure/layer_input/BiasAddBiasAdd3ResNet_entropy_closure/layer_input/MatMul:product:0AResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
0ResNet_entropy_closure/tf.math.softplus/SoftplusSoftplus3ResNet_entropy_closure/layer_input/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_0_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_0_layer_0/MatMulMatMul>ResNet_entropy_closure/tf.math.softplus/Softplus:activations:0DResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_0_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_0_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_0_layer_0/MatMul:product:0EResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
ResNet_entropy_closure/add/addAddV23ResNet_entropy_closure/layer_input/BiasAdd:output:07ResNet_entropy_closure/block_0_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_1/SoftplusSoftplus"ResNet_entropy_closure/add/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_1_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_1_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_1/Softplus:activations:0DResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_1_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_1_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_1_layer_0/MatMul:product:0EResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_1/addAddV2"ResNet_entropy_closure/add/add:z:07ResNet_entropy_closure/block_1_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_2/SoftplusSoftplus$ResNet_entropy_closure/add_1/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_2_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_2_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_2/Softplus:activations:0DResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_2_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_2_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_2_layer_0/MatMul:product:0EResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_2/addAddV2$ResNet_entropy_closure/add_1/add:z:07ResNet_entropy_closure/block_2_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_3/SoftplusSoftplus$ResNet_entropy_closure/add_2/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_3_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_3_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_3/Softplus:activations:0DResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_3_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_3_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_3_layer_0/MatMul:product:0EResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_3/addAddV2$ResNet_entropy_closure/add_2/add:z:07ResNet_entropy_closure/block_3_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_4/SoftplusSoftplus$ResNet_entropy_closure/add_3/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_4_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_4_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_4/Softplus:activations:0DResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_4_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_4_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_4_layer_0/MatMul:product:0EResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_4/addAddV2$ResNet_entropy_closure/add_3/add:z:07ResNet_entropy_closure/block_4_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_5/SoftplusSoftplus$ResNet_entropy_closure/add_4/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_5_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_5_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_5/Softplus:activations:0DResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_5_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_5_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_5_layer_0/MatMul:product:0EResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_5/addAddV2$ResNet_entropy_closure/add_4/add:z:07ResNet_entropy_closure/block_5_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
9ResNet_entropy_closure/dense_output/MatMul/ReadVariableOpReadVariableOpBresnet_entropy_closure_dense_output_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
*ResNet_entropy_closure/dense_output/MatMulMatMul$ResNet_entropy_closure/add_5/add:z:0AResNet_entropy_closure/dense_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
:ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOpReadVariableOpCresnet_entropy_closure_dense_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
+ResNet_entropy_closure/dense_output/BiasAddBiasAdd4ResNet_entropy_closure/dense_output/MatMul:product:0BResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������s
ones_like/ShapeShape4ResNet_entropy_closure/dense_output/BiasAdd:output:0*
T0*
_output_shapes
:T
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?w
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:����������
Egradient_tape/ResNet_entropy_closure/dense_output/BiasAdd/BiasAddGradBiasAddGradones_like:output:0*
T0*
_output_shapes
:�
?gradient_tape/ResNet_entropy_closure/dense_output/MatMul/MatMulMatMulones_like:output:0AResNet_entropy_closure/dense_output/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Agradient_tape/ResNet_entropy_closure/dense_output/MatMul/MatMul_1MatMul$ResNet_entropy_closure/add_5/add:z:0ones_like:output:0*
T0*
_output_shapes
:	�*
transpose_a(�
4gradient_tape/ResNet_entropy_closure/add_5/add/ShapeShape$ResNet_entropy_closure/add_4/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_5/add/Shape_1Shape7ResNet_entropy_closure/block_5_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_5/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_5/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_5/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_5/add/SumSumIgradient_tape/ResNet_entropy_closure/dense_output/MatMul/MatMul:product:0Igradient_tape/ResNet_entropy_closure/add_5/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_5/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_5/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_5/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_5/add/Sum_1SumIgradient_tape/ResNet_entropy_closure/dense_output/MatMul/MatMul:product:0Igradient_tape/ResNet_entropy_closure/add_5/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_5/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_5/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_5/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_5_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_5/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_5_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_5/add/Reshape_1:output:0DResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_5_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_5/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_5/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_5/SigmoidSigmoid$ResNet_entropy_closure/add_4/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_5/mulMulLgradient_tape/ResNet_entropy_closure/block_5_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_5/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddNAddN?gradient_tape/ResNet_entropy_closure/add_5/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_5/mul:z:0*
N*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_4/add/ShapeShape$ResNet_entropy_closure/add_3/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_4/add/Shape_1Shape7ResNet_entropy_closure/block_4_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_4/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_4/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_4/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_4/add/SumSum
AddN:sum:0Igradient_tape/ResNet_entropy_closure/add_4/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_4/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_4/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_4/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_4/add/Sum_1Sum
AddN:sum:0Igradient_tape/ResNet_entropy_closure/add_4/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_4/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_4/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_4/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_4_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_4/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_4_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_4/add/Reshape_1:output:0DResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_4_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_4/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_4/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_4/SigmoidSigmoid$ResNet_entropy_closure/add_3/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_4/mulMulLgradient_tape/ResNet_entropy_closure/block_4_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_4/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_1AddN?gradient_tape/ResNet_entropy_closure/add_4/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_4/mul:z:0*
N*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_3/add/ShapeShape$ResNet_entropy_closure/add_2/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_3/add/Shape_1Shape7ResNet_entropy_closure/block_3_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_3/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_3/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_3/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_3/add/SumSumAddN_1:sum:0Igradient_tape/ResNet_entropy_closure/add_3/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_3/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_3/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_3/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_3/add/Sum_1SumAddN_1:sum:0Igradient_tape/ResNet_entropy_closure/add_3/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_3/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_3/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_3/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_3_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_3/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_3_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_3/add/Reshape_1:output:0DResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_3_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_3/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_3/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_3/SigmoidSigmoid$ResNet_entropy_closure/add_2/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_3/mulMulLgradient_tape/ResNet_entropy_closure/block_3_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_3/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_2AddN?gradient_tape/ResNet_entropy_closure/add_3/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_3/mul:z:0*
N*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_2/add/ShapeShape$ResNet_entropy_closure/add_1/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_2/add/Shape_1Shape7ResNet_entropy_closure/block_2_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_2/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_2/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_2/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_2/add/SumSumAddN_2:sum:0Igradient_tape/ResNet_entropy_closure/add_2/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_2/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_2/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_2/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_2/add/Sum_1SumAddN_2:sum:0Igradient_tape/ResNet_entropy_closure/add_2/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_2/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_2/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_2/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_2_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_2/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_2_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_2/add/Reshape_1:output:0DResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_2_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_2/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_2/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_2/SigmoidSigmoid$ResNet_entropy_closure/add_1/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_2/mulMulLgradient_tape/ResNet_entropy_closure/block_2_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_2/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_3AddN?gradient_tape/ResNet_entropy_closure/add_2/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_2/mul:z:0*
N*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_1/add/ShapeShape"ResNet_entropy_closure/add/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_1/add/Shape_1Shape7ResNet_entropy_closure/block_1_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_1/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_1/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_1/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_1/add/SumSumAddN_3:sum:0Igradient_tape/ResNet_entropy_closure/add_1/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_1/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_1/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_1/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_1/add/Sum_1SumAddN_3:sum:0Igradient_tape/ResNet_entropy_closure/add_1/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_1/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_1/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_1/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_1_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_1/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_1_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_1/add/Reshape_1:output:0DResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_1_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_1/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_1/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_1/SigmoidSigmoid"ResNet_entropy_closure/add/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_1/mulMulLgradient_tape/ResNet_entropy_closure/block_1_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_1/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_4AddN?gradient_tape/ResNet_entropy_closure/add_1/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_1/mul:z:0*
N*
T0*(
_output_shapes
:�����������
2gradient_tape/ResNet_entropy_closure/add/add/ShapeShape3ResNet_entropy_closure/layer_input/BiasAdd:output:0*
T0*
_output_shapes
:�
4gradient_tape/ResNet_entropy_closure/add/add/Shape_1Shape7ResNet_entropy_closure/block_0_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Bgradient_tape/ResNet_entropy_closure/add/add/BroadcastGradientArgsBroadcastGradientArgs;gradient_tape/ResNet_entropy_closure/add/add/Shape:output:0=gradient_tape/ResNet_entropy_closure/add/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
0gradient_tape/ResNet_entropy_closure/add/add/SumSumAddN_4:sum:0Ggradient_tape/ResNet_entropy_closure/add/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
4gradient_tape/ResNet_entropy_closure/add/add/ReshapeReshape9gradient_tape/ResNet_entropy_closure/add/add/Sum:output:0;gradient_tape/ResNet_entropy_closure/add/add/Shape:output:0*
T0*(
_output_shapes
:�����������
2gradient_tape/ResNet_entropy_closure/add/add/Sum_1SumAddN_4:sum:0Ggradient_tape/ResNet_entropy_closure/add/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add/add/Reshape_1Reshape;gradient_tape/ResNet_entropy_closure/add/add/Sum_1:output:0=gradient_tape/ResNet_entropy_closure/add/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_0_layer_0/BiasAdd/BiasAddGradBiasAddGrad?gradient_tape/ResNet_entropy_closure/add/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_0_layer_0/MatMul/MatMulMatMul?gradient_tape/ResNet_entropy_closure/add/add/Reshape_1:output:0DResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_0_layer_0/MatMul/MatMul_1MatMul>ResNet_entropy_closure/tf.math.softplus/Softplus:activations:0?gradient_tape/ResNet_entropy_closure/add/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
=gradient_tape/ResNet_entropy_closure/tf.math.softplus/SigmoidSigmoid3ResNet_entropy_closure/layer_input/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
9gradient_tape/ResNet_entropy_closure/tf.math.softplus/mulMulLgradient_tape/ResNet_entropy_closure/block_0_layer_0/MatMul/MatMul:product:0Agradient_tape/ResNet_entropy_closure/tf.math.softplus/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_5AddN=gradient_tape/ResNet_entropy_closure/add/add/Reshape:output:0=gradient_tape/ResNet_entropy_closure/tf.math.softplus/mul:z:0*
N*
T0*(
_output_shapes
:�����������
Dgradient_tape/ResNet_entropy_closure/layer_input/BiasAdd/BiasAddGradBiasAddGradAddN_5:sum:0*
T0*
_output_shapes	
:��
>gradient_tape/ResNet_entropy_closure/layer_input/MatMul/MatMulMatMulAddN_5:sum:0@ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	*
transpose_b(�
@gradient_tape/ResNet_entropy_closure/layer_input/MatMul/MatMul_1MatMul;ResNet_entropy_closure/decorrelation_layer/MatMul:product:0AddN_5:sum:0*
T0*
_output_shapes
:		�*
transpose_a(�
Fgradient_tape/ResNet_entropy_closure/decorrelation_layer/MatMul/MatMulMatMulHgradient_tape/ResNet_entropy_closure/layer_input/MatMul/MatMul:product:0HResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	*
transpose_b(p
?gradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/ShapeShapex*
T0*
_output_shapes
:�
Agradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Shape_1ShapeBResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp:value:0*
T0*
_output_shapes
:�
Ogradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/BroadcastGradientArgsBroadcastGradientArgsHgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Shape:output:0Jgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Shape_1:output:0*2
_output_shapes 
:���������:����������
=gradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/SumSumPgradient_tape/ResNet_entropy_closure/decorrelation_layer/MatMul/MatMul:product:0Tgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
Agradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/ReshapeReshapeFgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Sum:output:0Hgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Shape:output:0*
T0*'
_output_shapes
:���������	�
CastCastJgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Reshape:output:0*

DstT0*

SrcT0*'
_output_shapes
:���������	�
checkedCheckNumericsCast:y:0*
T0*'
_output_shapes
:���������	*d
messageYWinput tensor checking error at alpha = Tensor("Cast:0", shape=(None, 9), dtype=float64)d
checkedandclipped/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped/MinimumMinimumchecked:output:0$checkedandclipped/Minimum/y:output:0*
T0*'
_output_shapes
:���������	\
checkedandclipped/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclippedMaximumcheckedandclipped/Minimum:z:0checkedandclipped/y:output:0*
T0*'
_output_shapes
:���������	d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSliceunknownstrided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:		�*

begin_mask*
end_maskX
Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:X
Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB: T
Tensordot/ShapeShapecheckedandclipped:z:0*
T0*
_output_shapes
:Y
Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2GatherV2Tensordot/Shape:output:0Tensordot/free:output:0 Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2_1GatherV2Tensordot/Shape:output:0Tensordot/axes:output:0"Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:Y
Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: n
Tensordot/ProdProdTensordot/GatherV2:output:0Tensordot/Const:output:0*
T0*
_output_shapes
: [
Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: t
Tensordot/Prod_1ProdTensordot/GatherV2_1:output:0Tensordot/Const_1:output:0*
T0*
_output_shapes
: W
Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concatConcatV2Tensordot/free:output:0Tensordot/axes:output:0Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:y
Tensordot/stackPackTensordot/Prod:output:0Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot/transpose	Transposecheckedandclipped:z:0Tensordot/concat:output:0*
T0*'
_output_shapes
:���������	�
Tensordot/ReshapeReshapeTensordot/transpose:y:0Tensordot/stack:output:0*
T0*0
_output_shapes
:�������������������
Tensordot/MatMulMatMulTensordot/Reshape:output:0strided_slice:output:0*
T0*(
_output_shapes
:����������\
Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�Y
Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concat_1ConcatV2Tensordot/GatherV2:output:0Tensordot/Const_2:output:0 Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
	TensordotReshapeTensordot/MatMul:product:0Tensordot/concat_1:output:0*
T0*(
_output_shapes
:����������Q
ExpExpTensordot:output:0*
T0*(
_output_shapes
:����������Z
Tensordot_1/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_1/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_1/ShapeShapeExp:y:0*
T0*
_output_shapes
:[
Tensordot_1/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2GatherV2Tensordot_1/Shape:output:0Tensordot_1/free:output:0"Tensordot_1/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_1/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2_1GatherV2Tensordot_1/Shape:output:0Tensordot_1/axes:output:0$Tensordot_1/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_1/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_1/ProdProdTensordot_1/GatherV2:output:0Tensordot_1/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_1/Prod_1ProdTensordot_1/GatherV2_1:output:0Tensordot_1/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concatConcatV2Tensordot_1/free:output:0Tensordot_1/axes:output:0 Tensordot_1/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_1/stackPackTensordot_1/Prod:output:0Tensordot_1/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_1/transpose	TransposeExp:y:0Tensordot_1/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_1/ReshapeReshapeTensordot_1/transpose:y:0Tensordot_1/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_1/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_1/transpose_1	Transposetensordot_1_b%Tensordot_1/transpose_1/perm:output:0*
T0*
_output_shapes
:	��
Tensordot_1/MatMulMatMulTensordot_1/Reshape:output:0Tensordot_1/transpose_1:y:0*
T0*'
_output_shapes
:���������]
Tensordot_1/Const_2Const*
_output_shapes
:*
dtype0*
valueB:[
Tensordot_1/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concat_1ConcatV2Tensordot_1/GatherV2:output:0Tensordot_1/Const_2:output:0"Tensordot_1/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_1ReshapeTensordot_1/MatMul:product:0Tensordot_1/concat_1:output:0*
T0*'
_output_shapes
:���������R
LogLogTensordot_1:output:0*
T0*'
_output_shapes
:���������E
NegNegLog:y:0*
T0*'
_output_shapes
:���������M
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :v
concatConcatV2Neg:y:0Cast:y:0concat/axis:output:0*
N*
T0*'
_output_shapes
:���������
�
	checked_1CheckNumericsconcat:output:0^checked*
T0*'
_output_shapes
:���������
*(
messageinput tensor checking errorf
checkedandclipped_1/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped_1/MinimumMinimumchecked_1:output:0&checkedandclipped_1/Minimum/y:output:0*
T0*'
_output_shapes
:���������
^
checkedandclipped_1/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclipped_1Maximumcheckedandclipped_1/Minimum:z:0checkedandclipped_1/y:output:0*
T0*'
_output_shapes
:���������
Z
Tensordot_2/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_2/freeConst*
_output_shapes
:*
dtype0*
valueB: X
Tensordot_2/ShapeShapecheckedandclipped_1:z:0*
T0*
_output_shapes
:[
Tensordot_2/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2GatherV2Tensordot_2/Shape:output:0Tensordot_2/free:output:0"Tensordot_2/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_2/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2_1GatherV2Tensordot_2/Shape:output:0Tensordot_2/axes:output:0$Tensordot_2/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_2/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_2/ProdProdTensordot_2/GatherV2:output:0Tensordot_2/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_2/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_2/Prod_1ProdTensordot_2/GatherV2_1:output:0Tensordot_2/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_2/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concatConcatV2Tensordot_2/free:output:0Tensordot_2/axes:output:0 Tensordot_2/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_2/stackPackTensordot_2/Prod:output:0Tensordot_2/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2/transpose	Transposecheckedandclipped_1:z:0Tensordot_2/concat:output:0*
T0*'
_output_shapes
:���������
�
Tensordot_2/ReshapeReshapeTensordot_2/transpose:y:0Tensordot_2/stack:output:0*
T0*0
_output_shapes
:������������������v
Tensordot_2/MatMulMatMulTensordot_2/Reshape:output:0unknown*
T0*(
_output_shapes
:����������^
Tensordot_2/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�[
Tensordot_2/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concat_1ConcatV2Tensordot_2/GatherV2:output:0Tensordot_2/Const_2:output:0"Tensordot_2/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2ReshapeTensordot_2/MatMul:product:0Tensordot_2/concat_1:output:0*
T0*(
_output_shapes
:����������U
Exp_1ExpTensordot_2:output:0*
T0*(
_output_shapes
:����������W
MulMul	Exp_1:y:0tensordot_1_b*
T0*(
_output_shapes
:����������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSliceunknownstrided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
:	
�*

begin_mask*
end_maskZ
Tensordot_3/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_3/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_3/ShapeShapeMul:z:0*
T0*
_output_shapes
:[
Tensordot_3/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2GatherV2Tensordot_3/Shape:output:0Tensordot_3/free:output:0"Tensordot_3/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_3/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2_1GatherV2Tensordot_3/Shape:output:0Tensordot_3/axes:output:0$Tensordot_3/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_3/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_3/ProdProdTensordot_3/GatherV2:output:0Tensordot_3/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_3/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_3/Prod_1ProdTensordot_3/GatherV2_1:output:0Tensordot_3/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_3/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concatConcatV2Tensordot_3/free:output:0Tensordot_3/axes:output:0 Tensordot_3/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_3/stackPackTensordot_3/Prod:output:0Tensordot_3/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_3/transpose	TransposeMul:z:0Tensordot_3/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_3/ReshapeReshapeTensordot_3/transpose:y:0Tensordot_3/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_3/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_3/transpose_1	Transposestrided_slice_1:output:0%Tensordot_3/transpose_1/perm:output:0*
T0*
_output_shapes
:	�
�
Tensordot_3/MatMulMatMulTensordot_3/Reshape:output:0Tensordot_3/transpose_1:y:0*
T0*'
_output_shapes
:���������
]
Tensordot_3/Const_2Const*
_output_shapes
:*
dtype0*
valueB:
[
Tensordot_3/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concat_1ConcatV2Tensordot_3/GatherV2:output:0Tensordot_3/Const_2:output:0"Tensordot_3/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_3ReshapeTensordot_3/MatMul:product:0Tensordot_3/concat_1:output:0*
T0*'
_output_shapes
:���������
X
Mul_1Mulmul_1_xconcat:output:0*
T0*'
_output_shapes
:���������
_
addAddV2Tensordot_3:output:0	Mul_1:z:0*
T0*'
_output_shapes
:���������
f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSliceadd:z:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������	*

begin_mask*
end_mask�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOpAresnet_entropy_closure_layer_input_matmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOpBresnet_entropy_closure_layer_input_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_0_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_0_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_1_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_1_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_2_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_2_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_3_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_3_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_4_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_4_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_5_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_5_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpBresnet_entropy_closure_dense_output_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
IdentityIdentity4ResNet_entropy_closure/dense_output/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������

Identity_1IdentityJgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Reshape:output:0^NoOp*
T0*'
_output_shapes
:���������	i

Identity_2Identitystrided_slice_2:output:0^NoOp*
T0*'
_output_shapes
:���������	�
NoOpNoOp>^ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOpA^ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp;^ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp:^ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp:^ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp9^ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp;^ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp^checked
^checked_16^dense_output/kernel/Regularizer/Square/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
2~
=ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp2�
@ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp@ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp2x
:ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp:ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp2v
9ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp9ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp2v
9ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp9ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp2t
8ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp8ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp2x
:ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp:ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2
checkedchecked2
	checked_1	checked_12n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp:J F
'
_output_shapes
:���������	

_user_specified_namex:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
�
A__forward_add_4_layer_call_and_return_conditional_losses_93976535
inputs_0

inputs_1_0
identity

inputs
inputs_1U
addAddV2inputs_0
inputs_1_0*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"
inputsinputs_0"
inputs_1
inputs_1_0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������*s
backward_function_nameYW__inference___backward_add_4_layer_call_and_return_conditional_losses_93976519_93976536:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_9_93980492N
?block_3_layer_0_bias_regularizer_square_readvariableop_resource:	�
identity��6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp?block_3_layer_0_bias_regularizer_square_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: f
IdentityIdentity(block_3_layer_0/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 
NoOpNoOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp
�
�
I__inference_layer_input_layer_call_and_return_conditional_losses_93980021

inputs1
matmul_readvariableop_resource:		�.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:		�*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������	: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
I__inference_layer_input_layer_call_and_return_conditional_losses_93975301

inputs1
matmul_readvariableop_resource:		�.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:		�*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������	: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
Q__inference_decorrelation_layer_layer_call_and_return_conditional_losses_93975275

inputs0
matmul_readvariableop_resource:		
identity��MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:		*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	_
IdentityIdentityMatMul:product:0^NoOp*
T0*'
_output_shapes
:���������	^
NoOpNoOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*(
_input_shapes
:���������	: 2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
N__inference_mean_shift_layer_layer_call_and_return_conditional_losses_93975264

inputs)
sub_readvariableop_resource:	
identity��sub/ReadVariableOpj
sub/ReadVariableOpReadVariableOpsub_readvariableop_resource*
_output_shapes
:	*
dtype0`
subSubinputssub/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	V
IdentityIdentitysub:z:0^NoOp*
T0*'
_output_shapes
:���������	[
NoOpNoOp^sub/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*(
_input_shapes
:���������	: 2(
sub/ReadVariableOpsub/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
G__forward_layer_input_layer_call_and_return_conditional_losses_93976776
inputs_01
matmul_readvariableop_resource:		�.
biasadd_readvariableop_resource:	�
identity
matmul_readvariableop

inputs��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:		�*
dtype0l
MatMulMatMulinputs_0MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"6
matmul_readvariableopMatMul/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������	: : *y
backward_function_name_]__inference___backward_layer_input_layer_call_and_return_conditional_losses_93976764_9397677720
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
T
(__inference_add_3_layer_call_fn_93980235
inputs_0
inputs_1
identity�
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_3_layer_call_and_return_conditional_losses_93975453a
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
�
�
J__inference_dense_output_layer_call_and_return_conditional_losses_93980382

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�5dense_output/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: _
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp6^dense_output/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976246
input_1'
mean_shift_layer_93976097:	.
decorrelation_layer_93976100:		'
layer_input_93976103:		�#
layer_input_93976105:	�,
block_0_layer_0_93976109:
��'
block_0_layer_0_93976111:	�,
block_1_layer_0_93976116:
��'
block_1_layer_0_93976118:	�,
block_2_layer_0_93976123:
��'
block_2_layer_0_93976125:	�,
block_3_layer_0_93976130:
��'
block_3_layer_0_93976132:	�,
block_4_layer_0_93976137:
��'
block_4_layer_0_93976139:	�,
block_5_layer_0_93976144:
��'
block_5_layer_0_93976146:	�(
dense_output_93976150:	�#
dense_output_93976152:
identity��'block_0_layer_0/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_1_layer_0/StatefulPartitionedCall�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_2_layer_0/StatefulPartitionedCall�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_3_layer_0/StatefulPartitionedCall�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_4_layer_0/StatefulPartitionedCall�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_5_layer_0/StatefulPartitionedCall�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�+decorrelation_layer/StatefulPartitionedCall�$dense_output/StatefulPartitionedCall�5dense_output/kernel/Regularizer/Square/ReadVariableOp�#layer_input/StatefulPartitionedCall�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�(mean_shift_layer/StatefulPartitionedCall�
(mean_shift_layer/StatefulPartitionedCallStatefulPartitionedCallinput_1mean_shift_layer_93976097*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *W
fRRP
N__inference_mean_shift_layer_layer_call_and_return_conditional_losses_93975264�
+decorrelation_layer/StatefulPartitionedCallStatefulPartitionedCall1mean_shift_layer/StatefulPartitionedCall:output:0decorrelation_layer_93976100*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *Z
fURS
Q__inference_decorrelation_layer_layer_call_and_return_conditional_losses_93975275�
#layer_input/StatefulPartitionedCallStatefulPartitionedCall4decorrelation_layer/StatefulPartitionedCall:output:0layer_input_93976103layer_input_93976105*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *R
fMRK
I__inference_layer_input_layer_call_and_return_conditional_losses_93975301�
tf.math.softplus/SoftplusSoftplus,layer_input/StatefulPartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_0_layer_0/StatefulPartitionedCallStatefulPartitionedCall'tf.math.softplus/Softplus:activations:0block_0_layer_0_93976109block_0_layer_0_93976111*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_0_layer_0_layer_call_and_return_conditional_losses_93975330�
add/PartitionedCallPartitionedCall,layer_input/StatefulPartitionedCall:output:00block_0_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__inference_add_layer_call_and_return_conditional_losses_93975342x
tf.math.softplus_1/SoftplusSoftplusadd/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_1_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_1/Softplus:activations:0block_1_layer_0_93976116block_1_layer_0_93976118*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_1_layer_0_layer_call_and_return_conditional_losses_93975367�
add_1/PartitionedCallPartitionedCalladd/PartitionedCall:output:00block_1_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_1_layer_call_and_return_conditional_losses_93975379z
tf.math.softplus_2/SoftplusSoftplusadd_1/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_2_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_2/Softplus:activations:0block_2_layer_0_93976123block_2_layer_0_93976125*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_2_layer_0_layer_call_and_return_conditional_losses_93975404�
add_2/PartitionedCallPartitionedCalladd_1/PartitionedCall:output:00block_2_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_2_layer_call_and_return_conditional_losses_93975416z
tf.math.softplus_3/SoftplusSoftplusadd_2/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_3_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_3/Softplus:activations:0block_3_layer_0_93976130block_3_layer_0_93976132*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_3_layer_0_layer_call_and_return_conditional_losses_93975441�
add_3/PartitionedCallPartitionedCalladd_2/PartitionedCall:output:00block_3_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_3_layer_call_and_return_conditional_losses_93975453z
tf.math.softplus_4/SoftplusSoftplusadd_3/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_4_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_4/Softplus:activations:0block_4_layer_0_93976137block_4_layer_0_93976139*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_4_layer_0_layer_call_and_return_conditional_losses_93975478�
add_4/PartitionedCallPartitionedCalladd_3/PartitionedCall:output:00block_4_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_4_layer_call_and_return_conditional_losses_93975490z
tf.math.softplus_5/SoftplusSoftplusadd_4/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_5_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_5/Softplus:activations:0block_5_layer_0_93976144block_5_layer_0_93976146*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_5_layer_0_layer_call_and_return_conditional_losses_93975515�
add_5/PartitionedCallPartitionedCalladd_4/PartitionedCall:output:00block_5_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_5_layer_call_and_return_conditional_losses_93975527�
$dense_output/StatefulPartitionedCallStatefulPartitionedCalladd_5/PartitionedCall:output:0dense_output_93976150dense_output_93976152*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *S
fNRL
J__inference_dense_output_layer_call_and_return_conditional_losses_93975545�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93976103*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93976105*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93976109* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93976111*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93976116* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93976118*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93976123* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93976125*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93976130* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93976132*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93976137* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93976139*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93976144* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93976146*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_output_93976150*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: |
IdentityIdentity-dense_output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������

NoOpNoOp(^block_0_layer_0/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_1_layer_0/StatefulPartitionedCall7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_2_layer_0/StatefulPartitionedCall7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_3_layer_0/StatefulPartitionedCall7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_4_layer_0/StatefulPartitionedCall7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_5_layer_0/StatefulPartitionedCall7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp,^decorrelation_layer/StatefulPartitionedCall%^dense_output/StatefulPartitionedCall6^dense_output/kernel/Regularizer/Square/ReadVariableOp$^layer_input/StatefulPartitionedCall3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp)^mean_shift_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 2R
'block_0_layer_0/StatefulPartitionedCall'block_0_layer_0/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_1_layer_0/StatefulPartitionedCall'block_1_layer_0/StatefulPartitionedCall2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_2_layer_0/StatefulPartitionedCall'block_2_layer_0/StatefulPartitionedCall2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_3_layer_0/StatefulPartitionedCall'block_3_layer_0/StatefulPartitionedCall2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_4_layer_0/StatefulPartitionedCall'block_4_layer_0/StatefulPartitionedCall2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_5_layer_0/StatefulPartitionedCall'block_5_layer_0/StatefulPartitionedCall2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2Z
+decorrelation_layer/StatefulPartitionedCall+decorrelation_layer/StatefulPartitionedCall2L
$dense_output/StatefulPartitionedCall$dense_output/StatefulPartitionedCall2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2J
#layer_input/StatefulPartitionedCall#layer_input/StatefulPartitionedCall2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp2T
(mean_shift_layer/StatefulPartitionedCall(mean_shift_layer/StatefulPartitionedCall:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1
�
�
M__inference_block_3_layer_0_layer_call_and_return_conditional_losses_93975441

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
K__forward_block_4_layer_0_layer_call_and_return_conditional_losses_93976555
inputs_02
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity
matmul_readvariableop

inputs��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0l
MatMulMatMulinputs_0MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"6
matmul_readvariableopMatMul/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : *}
backward_function_nameca__inference___backward_block_4_layer_0_layer_call_and_return_conditional_losses_93976543_9397655620
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
.__inference_layer_input_layer_call_fn_93979999

inputs
unknown:		�
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *R
fMRK
I__inference_layer_input_layer_call_and_return_conditional_losses_93975301p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������	: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
M__inference_block_0_layer_0_layer_call_and_return_conditional_losses_93980064

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
0__inference_sobolev_model_layer_call_fn_93977797
input_1
unknown:	
	unknown_0:		
	unknown_1:		�
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:
��

unknown_12:	�

unknown_13:
��

unknown_14:	�

unknown_15:	�

unknown_16:

unknown_17

unknown_18

unknown_19
identity

identity_1

identity_2��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19*!
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������	:���������	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93977697o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������	q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������	`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
R
&__inference_add_layer_call_fn_93980070
inputs_0
inputs_1
identity�
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__inference_add_layer_call_and_return_conditional_losses_93975342a
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
�
�
H__forward_dense_output_layer_call_and_return_conditional_losses_93976461
inputs_01
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity
matmul_readvariableop

inputs��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�5dense_output/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0k
MatMulMatMulinputs_0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: _
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp6^dense_output/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"6
matmul_readvariableopMatMul/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : *z
backward_function_name`^__inference___backward_dense_output_layer_call_and_return_conditional_losses_93976449_9397646220
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
T
(__inference_add_2_layer_call_fn_93980180
inputs_0
inputs_1
identity�
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_2_layer_call_and_return_conditional_losses_93975416a
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
�
�
M__inference_block_2_layer_0_layer_call_and_return_conditional_losses_93980174

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
W__inference___backward_add_2_layer_call_and_return_conditional_losses_93976617_93976634
placeholder#
gradients_add_grad_shape_inputs'
#gradients_add_grad_shape_1_inputs_1
identity

identity_1_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������g
gradients/add_grad/ShapeShapegradients_add_grad_shape_inputs*
T0*
_output_shapes
:m
gradients/add_grad/Shape_1Shape#gradients_add_grad_shape_1_inputs_1*
T0*
_output_shapes
:�
(gradients/add_grad/BroadcastGradientArgsBroadcastGradientArgs!gradients/add_grad/Shape:output:0#gradients/add_grad/Shape_1:output:0*2
_output_shapes 
:���������:����������
gradients/add_grad/SumSumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
gradients/add_grad/ReshapeReshapegradients/add_grad/Sum:output:0!gradients/add_grad/Shape:output:0*
T0*(
_output_shapes
:�����������
gradients/add_grad/Sum_1Sumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
gradients/add_grad/Reshape_1Reshape!gradients/add_grad/Sum_1:output:0#gradients/add_grad/Shape_1:output:0*
T0*(
_output_shapes
:����������l
IdentityIdentity#gradients/add_grad/Reshape:output:0*
T0*(
_output_shapes
:����������p

Identity_1Identity%gradients/add_grad/Reshape_1:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������:����������:����������*\
forward_function_nameCA__forward_add_2_layer_call_and_return_conditional_losses_93976633:. *
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������
�
k
A__inference_add_layer_call_and_return_conditional_losses_93975342

inputs
inputs_1
identityQ
addAddV2inputsinputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
^__inference___backward_dense_output_layer_call_and_return_conditional_losses_93976449_93976462
placeholder6
2gradients_matmul_grad_matmul_matmul_readvariableop)
%gradients_matmul_grad_matmul_1_inputs
identity

identity_1

identity_2^
gradients/grad_ys_0Identityplaceholder*
T0*'
_output_shapes
:���������t
"gradients/BiasAdd_grad/BiasAddGradBiasAddGradgradients/grad_ys_0:output:0*
T0*
_output_shapes
:�
gradients/MatMul_grad/MatMulMatMulgradients/grad_ys_0:output:02gradients_matmul_grad_matmul_matmul_readvariableop*
T0*(
_output_shapes
:����������*
transpose_b(�
gradients/MatMul_grad/MatMul_1MatMul%gradients_matmul_grad_matmul_1_inputsgradients/grad_ys_0:output:0*
T0*
_output_shapes
:	�*
transpose_a(o
IdentityIdentity&gradients/MatMul_grad/MatMul:product:0*
T0*(
_output_shapes
:����������j

Identity_1Identity(gradients/MatMul_grad/MatMul_1:product:0*
T0*
_output_shapes
:	�h

Identity_2Identity+gradients/BiasAdd_grad/BiasAddGrad:output:0*
T0*
_output_shapes
:"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*E
_input_shapes4
2:���������:	�:����������*c
forward_function_nameJH__forward_dense_output_layer_call_and_return_conditional_losses_93976461:- )
'
_output_shapes
:���������:%!

_output_shapes
:	�:.*
(
_output_shapes
:����������
�
�
a__inference___backward_block_5_layer_0_layer_call_and_return_conditional_losses_93976494_93976507
placeholder6
2gradients_matmul_grad_matmul_matmul_readvariableop)
%gradients_matmul_grad_matmul_1_inputs
identity

identity_1

identity_2_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������u
"gradients/BiasAdd_grad/BiasAddGradBiasAddGradgradients/grad_ys_0:output:0*
T0*
_output_shapes	
:��
gradients/MatMul_grad/MatMulMatMulgradients/grad_ys_0:output:02gradients_matmul_grad_matmul_matmul_readvariableop*
T0*(
_output_shapes
:����������*
transpose_b(�
gradients/MatMul_grad/MatMul_1MatMul%gradients_matmul_grad_matmul_1_inputsgradients/grad_ys_0:output:0*
T0* 
_output_shapes
:
��*
transpose_a(o
IdentityIdentity&gradients/MatMul_grad/MatMul:product:0*
T0*(
_output_shapes
:����������k

Identity_1Identity(gradients/MatMul_grad/MatMul_1:product:0*
T0* 
_output_shapes
:
��i

Identity_2Identity+gradients/BiasAdd_grad/BiasAddGrad:output:0*
T0*
_output_shapes	
:�"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*G
_input_shapes6
4:����������:
��:����������*f
forward_function_nameMK__forward_block_5_layer_0_layer_call_and_return_conditional_losses_93976506:. *
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������
�
�
0__inference_sobolev_model_layer_call_fn_93978684
x
unknown:	
	unknown_0:		
	unknown_1:		�
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:
��

unknown_12:	�

unknown_13:
��

unknown_14:	�

unknown_15:	�

unknown_16:

unknown_17

unknown_18

unknown_19
identity

identity_1

identity_2��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallxunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19*!
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������	:���������	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93977697o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������	q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������	`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
22
StatefulPartitionedCallStatefulPartitionedCall:J F
'
_output_shapes
:���������	

_user_specified_namex:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
T
(__inference_add_5_layer_call_fn_93980345
inputs_0
inputs_1
identity�
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_5_layer_call_and_return_conditional_losses_93975527a
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
ۆ
�
!__inference__traced_save_93980800
file_prefix(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop'
#savev2_variable_read_readvariableop)
%savev2_variable_1_read_readvariableop1
-savev2_layer_input_kernel_read_readvariableop/
+savev2_layer_input_bias_read_readvariableop5
1savev2_block_0_layer_0_kernel_read_readvariableop3
/savev2_block_0_layer_0_bias_read_readvariableop5
1savev2_block_1_layer_0_kernel_read_readvariableop3
/savev2_block_1_layer_0_bias_read_readvariableop5
1savev2_block_2_layer_0_kernel_read_readvariableop3
/savev2_block_2_layer_0_bias_read_readvariableop5
1savev2_block_3_layer_0_kernel_read_readvariableop3
/savev2_block_3_layer_0_bias_read_readvariableop5
1savev2_block_4_layer_0_kernel_read_readvariableop3
/savev2_block_4_layer_0_bias_read_readvariableop5
1savev2_block_5_layer_0_kernel_read_readvariableop3
/savev2_block_5_layer_0_bias_read_readvariableop2
.savev2_dense_output_kernel_read_readvariableop0
,savev2_dense_output_bias_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop&
"savev2_total_2_read_readvariableop&
"savev2_count_2_read_readvariableop&
"savev2_total_3_read_readvariableop&
"savev2_count_3_read_readvariableop&
"savev2_total_4_read_readvariableop&
"savev2_count_4_read_readvariableop&
"savev2_total_5_read_readvariableop&
"savev2_count_5_read_readvariableop&
"savev2_total_6_read_readvariableop&
"savev2_count_6_read_readvariableop&
"savev2_total_7_read_readvariableop&
"savev2_count_7_read_readvariableop&
"savev2_total_8_read_readvariableop&
"savev2_count_8_read_readvariableop&
"savev2_total_9_read_readvariableop&
"savev2_count_9_read_readvariableop8
4savev2_adam_layer_input_kernel_m_read_readvariableop6
2savev2_adam_layer_input_bias_m_read_readvariableop<
8savev2_adam_block_0_layer_0_kernel_m_read_readvariableop:
6savev2_adam_block_0_layer_0_bias_m_read_readvariableop<
8savev2_adam_block_1_layer_0_kernel_m_read_readvariableop:
6savev2_adam_block_1_layer_0_bias_m_read_readvariableop<
8savev2_adam_block_2_layer_0_kernel_m_read_readvariableop:
6savev2_adam_block_2_layer_0_bias_m_read_readvariableop<
8savev2_adam_block_3_layer_0_kernel_m_read_readvariableop:
6savev2_adam_block_3_layer_0_bias_m_read_readvariableop<
8savev2_adam_block_4_layer_0_kernel_m_read_readvariableop:
6savev2_adam_block_4_layer_0_bias_m_read_readvariableop<
8savev2_adam_block_5_layer_0_kernel_m_read_readvariableop:
6savev2_adam_block_5_layer_0_bias_m_read_readvariableop9
5savev2_adam_dense_output_kernel_m_read_readvariableop7
3savev2_adam_dense_output_bias_m_read_readvariableop8
4savev2_adam_layer_input_kernel_v_read_readvariableop6
2savev2_adam_layer_input_bias_v_read_readvariableop<
8savev2_adam_block_0_layer_0_kernel_v_read_readvariableop:
6savev2_adam_block_0_layer_0_bias_v_read_readvariableop<
8savev2_adam_block_1_layer_0_kernel_v_read_readvariableop:
6savev2_adam_block_1_layer_0_bias_v_read_readvariableop<
8savev2_adam_block_2_layer_0_kernel_v_read_readvariableop:
6savev2_adam_block_2_layer_0_bias_v_read_readvariableop<
8savev2_adam_block_3_layer_0_kernel_v_read_readvariableop:
6savev2_adam_block_3_layer_0_bias_v_read_readvariableop<
8savev2_adam_block_4_layer_0_kernel_v_read_readvariableop:
6savev2_adam_block_4_layer_0_bias_v_read_readvariableop<
8savev2_adam_block_5_layer_0_kernel_v_read_readvariableop:
6savev2_adam_block_5_layer_0_bias_v_read_readvariableop9
5savev2_adam_dense_output_kernel_v_read_readvariableop7
3savev2_adam_dense_output_bias_v_read_readvariableop
savev2_const_3

identity_1��MergeV2Checkpointsw
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*Z
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.parta
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: f

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: L

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :f
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : �
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �"
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:L*
dtype0*�!
value�!B�!LB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB'variables/14/.ATTRIBUTES/VARIABLE_VALUEB'variables/15/.ATTRIBUTES/VARIABLE_VALUEB'variables/16/.ATTRIBUTES/VARIABLE_VALUEB'variables/17/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/3/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/3/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/4/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/4/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/5/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/5/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/6/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/6/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/7/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/7/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/8/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/8/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/9/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/9/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/16/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/17/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/16/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/17/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:L*
dtype0*�
value�B�LB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop#savev2_variable_read_readvariableop%savev2_variable_1_read_readvariableop-savev2_layer_input_kernel_read_readvariableop+savev2_layer_input_bias_read_readvariableop1savev2_block_0_layer_0_kernel_read_readvariableop/savev2_block_0_layer_0_bias_read_readvariableop1savev2_block_1_layer_0_kernel_read_readvariableop/savev2_block_1_layer_0_bias_read_readvariableop1savev2_block_2_layer_0_kernel_read_readvariableop/savev2_block_2_layer_0_bias_read_readvariableop1savev2_block_3_layer_0_kernel_read_readvariableop/savev2_block_3_layer_0_bias_read_readvariableop1savev2_block_4_layer_0_kernel_read_readvariableop/savev2_block_4_layer_0_bias_read_readvariableop1savev2_block_5_layer_0_kernel_read_readvariableop/savev2_block_5_layer_0_bias_read_readvariableop.savev2_dense_output_kernel_read_readvariableop,savev2_dense_output_bias_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop"savev2_total_2_read_readvariableop"savev2_count_2_read_readvariableop"savev2_total_3_read_readvariableop"savev2_count_3_read_readvariableop"savev2_total_4_read_readvariableop"savev2_count_4_read_readvariableop"savev2_total_5_read_readvariableop"savev2_count_5_read_readvariableop"savev2_total_6_read_readvariableop"savev2_count_6_read_readvariableop"savev2_total_7_read_readvariableop"savev2_count_7_read_readvariableop"savev2_total_8_read_readvariableop"savev2_count_8_read_readvariableop"savev2_total_9_read_readvariableop"savev2_count_9_read_readvariableop4savev2_adam_layer_input_kernel_m_read_readvariableop2savev2_adam_layer_input_bias_m_read_readvariableop8savev2_adam_block_0_layer_0_kernel_m_read_readvariableop6savev2_adam_block_0_layer_0_bias_m_read_readvariableop8savev2_adam_block_1_layer_0_kernel_m_read_readvariableop6savev2_adam_block_1_layer_0_bias_m_read_readvariableop8savev2_adam_block_2_layer_0_kernel_m_read_readvariableop6savev2_adam_block_2_layer_0_bias_m_read_readvariableop8savev2_adam_block_3_layer_0_kernel_m_read_readvariableop6savev2_adam_block_3_layer_0_bias_m_read_readvariableop8savev2_adam_block_4_layer_0_kernel_m_read_readvariableop6savev2_adam_block_4_layer_0_bias_m_read_readvariableop8savev2_adam_block_5_layer_0_kernel_m_read_readvariableop6savev2_adam_block_5_layer_0_bias_m_read_readvariableop5savev2_adam_dense_output_kernel_m_read_readvariableop3savev2_adam_dense_output_bias_m_read_readvariableop4savev2_adam_layer_input_kernel_v_read_readvariableop2savev2_adam_layer_input_bias_v_read_readvariableop8savev2_adam_block_0_layer_0_kernel_v_read_readvariableop6savev2_adam_block_0_layer_0_bias_v_read_readvariableop8savev2_adam_block_1_layer_0_kernel_v_read_readvariableop6savev2_adam_block_1_layer_0_bias_v_read_readvariableop8savev2_adam_block_2_layer_0_kernel_v_read_readvariableop6savev2_adam_block_2_layer_0_bias_v_read_readvariableop8savev2_adam_block_3_layer_0_kernel_v_read_readvariableop6savev2_adam_block_3_layer_0_bias_v_read_readvariableop8savev2_adam_block_4_layer_0_kernel_v_read_readvariableop6savev2_adam_block_4_layer_0_bias_v_read_readvariableop8savev2_adam_block_5_layer_0_kernel_v_read_readvariableop6savev2_adam_block_5_layer_0_bias_v_read_readvariableop5savev2_adam_dense_output_kernel_v_read_readvariableop3savev2_adam_dense_output_bias_v_read_readvariableopsavev2_const_3"/device:CPU:0*
_output_shapes
 *Z
dtypesP
N2L	�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 f
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: Q

Identity_1IdentityIdentity:output:0^NoOp*
T0*
_output_shapes
: [
NoOpNoOp^MergeV2Checkpoints*"
_acd_function_control_output(*
_output_shapes
 "!

identity_1Identity_1:output:0*�
_input_shapes�
�: : : : : : :	:		:		�:�:
��:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�:: : : : : : : : : : : : : : : : : : : : :		�:�:
��:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�::		�:�:
��:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�:: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: : 

_output_shapes
:	:$ 

_output_shapes

:		:%!

_output_shapes
:		�:!	

_output_shapes	
:�:&
"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:%!

_output_shapes
:	�: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: : 

_output_shapes
: :!

_output_shapes
: :"

_output_shapes
: :#

_output_shapes
: :$

_output_shapes
: :%

_output_shapes
: :&

_output_shapes
: :'

_output_shapes
: :(

_output_shapes
: :)

_output_shapes
: :*

_output_shapes
: :+

_output_shapes
: :%,!

_output_shapes
:		�:!-

_output_shapes	
:�:&."
 
_output_shapes
:
��:!/

_output_shapes	
:�:&0"
 
_output_shapes
:
��:!1

_output_shapes	
:�:&2"
 
_output_shapes
:
��:!3

_output_shapes	
:�:&4"
 
_output_shapes
:
��:!5

_output_shapes	
:�:&6"
 
_output_shapes
:
��:!7

_output_shapes	
:�:&8"
 
_output_shapes
:
��:!9

_output_shapes	
:�:%:!

_output_shapes
:	�: ;

_output_shapes
::%<!

_output_shapes
:		�:!=

_output_shapes	
:�:&>"
 
_output_shapes
:
��:!?

_output_shapes	
:�:&@"
 
_output_shapes
:
��:!A

_output_shapes	
:�:&B"
 
_output_shapes
:
��:!C

_output_shapes	
:�:&D"
 
_output_shapes
:
��:!E

_output_shapes	
:�:&F"
 
_output_shapes
:
��:!G

_output_shapes	
:�:&H"
 
_output_shapes
:
��:!I

_output_shapes	
:�:%J!

_output_shapes
:	�: K

_output_shapes
::L

_output_shapes
: 
��
�
R__forward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976847

inputs'
mean_shift_layer_93975265:	.
decorrelation_layer_93975276:		'
layer_input_93975302:		�#
layer_input_93975304:	�,
block_0_layer_0_93975331:
��'
block_0_layer_0_93975333:	�,
block_1_layer_0_93975368:
��'
block_1_layer_0_93975370:	�,
block_2_layer_0_93975405:
��'
block_2_layer_0_93975407:	�,
block_3_layer_0_93975442:
��'
block_3_layer_0_93975444:	�,
block_4_layer_0_93975479:
��'
block_4_layer_0_93975481:	�,
block_5_layer_0_93975516:
��'
block_5_layer_0_93975518:	�(
dense_output_93975546:	�#
dense_output_93975548:
identity(
$dense_output_statefulpartitionedcall*
&dense_output_statefulpartitionedcall_0
add_5_partitionedcall
add_5_partitionedcall_0+
'block_5_layer_0_statefulpartitionedcall-
)block_5_layer_0_statefulpartitionedcall_0
add_4_partitionedcall
add_4_partitionedcall_0
add_4_partitionedcall_1+
'block_4_layer_0_statefulpartitionedcall-
)block_4_layer_0_statefulpartitionedcall_0
add_3_partitionedcall
add_3_partitionedcall_0
add_3_partitionedcall_1+
'block_3_layer_0_statefulpartitionedcall-
)block_3_layer_0_statefulpartitionedcall_0
add_2_partitionedcall
add_2_partitionedcall_0
add_2_partitionedcall_1+
'block_2_layer_0_statefulpartitionedcall-
)block_2_layer_0_statefulpartitionedcall_0
add_1_partitionedcall
add_1_partitionedcall_0
add_1_partitionedcall_1+
'block_1_layer_0_statefulpartitionedcall-
)block_1_layer_0_statefulpartitionedcall_0
add_partitionedcall
add_partitionedcall_0
add_partitionedcall_1+
'block_0_layer_0_statefulpartitionedcall-
)block_0_layer_0_statefulpartitionedcall_0'
#layer_input_statefulpartitionedcall)
%layer_input_statefulpartitionedcall_0)
%layer_input_statefulpartitionedcall_1/
+decorrelation_layer_statefulpartitionedcall1
-decorrelation_layer_statefulpartitionedcall_0,
(mean_shift_layer_statefulpartitionedcall.
*mean_shift_layer_statefulpartitionedcall_0��'block_0_layer_0/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_1_layer_0/StatefulPartitionedCall�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_2_layer_0/StatefulPartitionedCall�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_3_layer_0/StatefulPartitionedCall�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_4_layer_0/StatefulPartitionedCall�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_5_layer_0/StatefulPartitionedCall�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�+decorrelation_layer/StatefulPartitionedCall�$dense_output/StatefulPartitionedCall�5dense_output/kernel/Regularizer/Square/ReadVariableOp�#layer_input/StatefulPartitionedCall�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�(mean_shift_layer/StatefulPartitionedCall�
(mean_shift_layer/StatefulPartitionedCallStatefulPartitionedCallinputsmean_shift_layer_93975265*
Tin
2*
Tout
2*
_collective_manager_ids
 *@
_output_shapes.
,:���������	:���������	:	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *U
fPRN
L__forward_mean_shift_layer_layer_call_and_return_conditional_losses_93976820�
+decorrelation_layer/StatefulPartitionedCallStatefulPartitionedCall1mean_shift_layer/StatefulPartitionedCall:output:0decorrelation_layer_93975276*
Tin
2*
Tout
2*
_collective_manager_ids
 *D
_output_shapes2
0:���������	:		:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *X
fSRQ
O__forward_decorrelation_layer_layer_call_and_return_conditional_losses_93976795�
#layer_input/StatefulPartitionedCallStatefulPartitionedCall4decorrelation_layer/StatefulPartitionedCall:output:0layer_input_93975302layer_input_93975304*
Tin
2*
Tout
2*
_collective_manager_ids
 *F
_output_shapes4
2:����������:		�:���������	*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *P
fKRI
G__forward_layer_input_layer_call_and_return_conditional_losses_93976776�
tf.math.softplus/SoftplusSoftplus,layer_input/StatefulPartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_0_layer_0/StatefulPartitionedCallStatefulPartitionedCall'tf.math.softplus/Softplus:activations:0block_0_layer_0_93975331block_0_layer_0_93975333*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_0_layer_0_layer_call_and_return_conditional_losses_93976751�
add/PartitionedCallPartitionedCall,layer_input/StatefulPartitionedCall:output:00block_0_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *H
fCRA
?__forward_add_layer_call_and_return_conditional_losses_93976731x
tf.math.softplus_1/SoftplusSoftplusadd/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_1_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_1/Softplus:activations:0block_1_layer_0_93975368block_1_layer_0_93975370*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_1_layer_0_layer_call_and_return_conditional_losses_93976702�
add_1/PartitionedCallPartitionedCalladd/PartitionedCall:output:00block_1_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_1_layer_call_and_return_conditional_losses_93976682z
tf.math.softplus_2/SoftplusSoftplusadd_1/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_2_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_2/Softplus:activations:0block_2_layer_0_93975405block_2_layer_0_93975407*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_2_layer_0_layer_call_and_return_conditional_losses_93976653�
add_2/PartitionedCallPartitionedCalladd_1/PartitionedCall:output:00block_2_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_2_layer_call_and_return_conditional_losses_93976633z
tf.math.softplus_3/SoftplusSoftplusadd_2/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_3_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_3/Softplus:activations:0block_3_layer_0_93975442block_3_layer_0_93975444*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_3_layer_0_layer_call_and_return_conditional_losses_93976604�
add_3/PartitionedCallPartitionedCalladd_2/PartitionedCall:output:00block_3_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_3_layer_call_and_return_conditional_losses_93976584z
tf.math.softplus_4/SoftplusSoftplusadd_3/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_4_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_4/Softplus:activations:0block_4_layer_0_93975479block_4_layer_0_93975481*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_4_layer_0_layer_call_and_return_conditional_losses_93976555�
add_4/PartitionedCallPartitionedCalladd_3/PartitionedCall:output:00block_4_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_4_layer_call_and_return_conditional_losses_93976535z
tf.math.softplus_5/SoftplusSoftplusadd_4/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_5_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_5/Softplus:activations:0block_5_layer_0_93975516block_5_layer_0_93975518*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_5_layer_0_layer_call_and_return_conditional_losses_93976506�
add_5/PartitionedCallPartitionedCalladd_4/PartitionedCall:output:00block_5_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_5_layer_call_and_return_conditional_losses_93976486�
$dense_output/StatefulPartitionedCallStatefulPartitionedCalladd_5/PartitionedCall:output:0dense_output_93975546dense_output_93975548*
Tin
2*
Tout
2*
_collective_manager_ids
 *F
_output_shapes4
2:���������:	�:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *Q
fLRJ
H__forward_dense_output_layer_call_and_return_conditional_losses_93976461�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93975302*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93975304*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93975331* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93975333*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93975368* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93975370*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93975405* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93975407*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93975442* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93975444*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93975479* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93975481*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93975516* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93975518*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_output_93975546*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: |
IdentityIdentity-dense_output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������

NoOpNoOp(^block_0_layer_0/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_1_layer_0/StatefulPartitionedCall7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_2_layer_0/StatefulPartitionedCall7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_3_layer_0/StatefulPartitionedCall7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_4_layer_0/StatefulPartitionedCall7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_5_layer_0/StatefulPartitionedCall7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp,^decorrelation_layer/StatefulPartitionedCall%^dense_output/StatefulPartitionedCall6^dense_output/kernel/Regularizer/Square/ReadVariableOp$^layer_input/StatefulPartitionedCall3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp)^mean_shift_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "7
add_1_partitionedcalladd_1/PartitionedCall:output:0"9
add_1_partitionedcall_0add_1/PartitionedCall:output:1"9
add_1_partitionedcall_1add_1/PartitionedCall:output:2"7
add_2_partitionedcalladd_2/PartitionedCall:output:0"9
add_2_partitionedcall_0add_2/PartitionedCall:output:1"9
add_2_partitionedcall_1add_2/PartitionedCall:output:2"7
add_3_partitionedcalladd_3/PartitionedCall:output:0"9
add_3_partitionedcall_0add_3/PartitionedCall:output:1"9
add_3_partitionedcall_1add_3/PartitionedCall:output:2"7
add_4_partitionedcalladd_4/PartitionedCall:output:0"9
add_4_partitionedcall_0add_4/PartitionedCall:output:1"9
add_4_partitionedcall_1add_4/PartitionedCall:output:2"7
add_5_partitionedcalladd_5/PartitionedCall:output:1"9
add_5_partitionedcall_0add_5/PartitionedCall:output:2"3
add_partitionedcalladd/PartitionedCall:output:0"5
add_partitionedcall_0add/PartitionedCall:output:1"5
add_partitionedcall_1add/PartitionedCall:output:2"[
'block_0_layer_0_statefulpartitionedcall0block_0_layer_0/StatefulPartitionedCall:output:1"]
)block_0_layer_0_statefulpartitionedcall_00block_0_layer_0/StatefulPartitionedCall:output:2"[
'block_1_layer_0_statefulpartitionedcall0block_1_layer_0/StatefulPartitionedCall:output:1"]
)block_1_layer_0_statefulpartitionedcall_00block_1_layer_0/StatefulPartitionedCall:output:2"[
'block_2_layer_0_statefulpartitionedcall0block_2_layer_0/StatefulPartitionedCall:output:1"]
)block_2_layer_0_statefulpartitionedcall_00block_2_layer_0/StatefulPartitionedCall:output:2"[
'block_3_layer_0_statefulpartitionedcall0block_3_layer_0/StatefulPartitionedCall:output:1"]
)block_3_layer_0_statefulpartitionedcall_00block_3_layer_0/StatefulPartitionedCall:output:2"[
'block_4_layer_0_statefulpartitionedcall0block_4_layer_0/StatefulPartitionedCall:output:1"]
)block_4_layer_0_statefulpartitionedcall_00block_4_layer_0/StatefulPartitionedCall:output:2"[
'block_5_layer_0_statefulpartitionedcall0block_5_layer_0/StatefulPartitionedCall:output:1"]
)block_5_layer_0_statefulpartitionedcall_00block_5_layer_0/StatefulPartitionedCall:output:2"c
+decorrelation_layer_statefulpartitionedcall4decorrelation_layer/StatefulPartitionedCall:output:1"e
-decorrelation_layer_statefulpartitionedcall_04decorrelation_layer/StatefulPartitionedCall:output:2"U
$dense_output_statefulpartitionedcall-dense_output/StatefulPartitionedCall:output:1"W
&dense_output_statefulpartitionedcall_0-dense_output/StatefulPartitionedCall:output:2"
identityIdentity:output:0"S
#layer_input_statefulpartitionedcall,layer_input/StatefulPartitionedCall:output:0"U
%layer_input_statefulpartitionedcall_0,layer_input/StatefulPartitionedCall:output:1"U
%layer_input_statefulpartitionedcall_1,layer_input/StatefulPartitionedCall:output:2"]
(mean_shift_layer_statefulpartitionedcall1mean_shift_layer/StatefulPartitionedCall:output:1"_
*mean_shift_layer_statefulpartitionedcall_01mean_shift_layer/StatefulPartitionedCall:output:2*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : *�
backward_function_namejh__inference___backward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976445_939768482R
'block_0_layer_0/StatefulPartitionedCall'block_0_layer_0/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_1_layer_0/StatefulPartitionedCall'block_1_layer_0/StatefulPartitionedCall2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_2_layer_0/StatefulPartitionedCall'block_2_layer_0/StatefulPartitionedCall2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_3_layer_0/StatefulPartitionedCall'block_3_layer_0/StatefulPartitionedCall2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_4_layer_0/StatefulPartitionedCall'block_4_layer_0/StatefulPartitionedCall2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_5_layer_0/StatefulPartitionedCall'block_5_layer_0/StatefulPartitionedCall2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2Z
+decorrelation_layer/StatefulPartitionedCall+decorrelation_layer/StatefulPartitionedCall2L
$dense_output/StatefulPartitionedCall$dense_output/StatefulPartitionedCall2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2J
#layer_input/StatefulPartitionedCall#layer_input/StatefulPartitionedCall2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp2T
(mean_shift_layer/StatefulPartitionedCall(mean_shift_layer/StatefulPartitionedCall:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
a__inference___backward_block_0_layer_0_layer_call_and_return_conditional_losses_93976739_93976752
placeholder6
2gradients_matmul_grad_matmul_matmul_readvariableop)
%gradients_matmul_grad_matmul_1_inputs
identity

identity_1

identity_2_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������u
"gradients/BiasAdd_grad/BiasAddGradBiasAddGradgradients/grad_ys_0:output:0*
T0*
_output_shapes	
:��
gradients/MatMul_grad/MatMulMatMulgradients/grad_ys_0:output:02gradients_matmul_grad_matmul_matmul_readvariableop*
T0*(
_output_shapes
:����������*
transpose_b(�
gradients/MatMul_grad/MatMul_1MatMul%gradients_matmul_grad_matmul_1_inputsgradients/grad_ys_0:output:0*
T0* 
_output_shapes
:
��*
transpose_a(o
IdentityIdentity&gradients/MatMul_grad/MatMul:product:0*
T0*(
_output_shapes
:����������k

Identity_1Identity(gradients/MatMul_grad/MatMul_1:product:0*
T0* 
_output_shapes
:
��i

Identity_2Identity+gradients/BiasAdd_grad/BiasAddGrad:output:0*
T0*
_output_shapes	
:�"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*G
_input_shapes6
4:����������:
��:����������*f
forward_function_nameMK__forward_block_0_layer_0_layer_call_and_return_conditional_losses_93976751:. *
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������
�
�
K__forward_block_3_layer_0_layer_call_and_return_conditional_losses_93976604
inputs_02
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity
matmul_readvariableop

inputs��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0l
MatMulMatMulinputs_0MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"6
matmul_readvariableopMatMul/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : *}
backward_function_nameca__inference___backward_block_3_layer_0_layer_call_and_return_conditional_losses_93976592_9397660520
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
M__inference_block_5_layer_0_layer_call_and_return_conditional_losses_93980339

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
M__inference_block_5_layer_0_layer_call_and_return_conditional_losses_93975515

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
W__inference___backward_add_3_layer_call_and_return_conditional_losses_93976568_93976585
placeholder#
gradients_add_grad_shape_inputs'
#gradients_add_grad_shape_1_inputs_1
identity

identity_1_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������g
gradients/add_grad/ShapeShapegradients_add_grad_shape_inputs*
T0*
_output_shapes
:m
gradients/add_grad/Shape_1Shape#gradients_add_grad_shape_1_inputs_1*
T0*
_output_shapes
:�
(gradients/add_grad/BroadcastGradientArgsBroadcastGradientArgs!gradients/add_grad/Shape:output:0#gradients/add_grad/Shape_1:output:0*2
_output_shapes 
:���������:����������
gradients/add_grad/SumSumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
gradients/add_grad/ReshapeReshapegradients/add_grad/Sum:output:0!gradients/add_grad/Shape:output:0*
T0*(
_output_shapes
:�����������
gradients/add_grad/Sum_1Sumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
gradients/add_grad/Reshape_1Reshape!gradients/add_grad/Sum_1:output:0#gradients/add_grad/Shape_1:output:0*
T0*(
_output_shapes
:����������l
IdentityIdentity#gradients/add_grad/Reshape:output:0*
T0*(
_output_shapes
:����������p

Identity_1Identity%gradients/add_grad/Reshape_1:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������:����������:����������*\
forward_function_nameCA__forward_add_3_layer_call_and_return_conditional_losses_93976584:. *
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������
�
�
__inference_loss_fn_3_93980426N
?block_0_layer_0_bias_regularizer_square_readvariableop_resource:	�
identity��6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp?block_0_layer_0_bias_regularizer_square_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: f
IdentityIdentity(block_0_layer_0/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 
NoOpNoOp7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp
�
�
2__inference_block_5_layer_0_layer_call_fn_93980317

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_5_layer_0_layer_call_and_return_conditional_losses_93975515p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93979950

inputs:
,mean_shift_layer_sub_readvariableop_resource:	D
2decorrelation_layer_matmul_readvariableop_resource:		=
*layer_input_matmul_readvariableop_resource:		�:
+layer_input_biasadd_readvariableop_resource:	�B
.block_0_layer_0_matmul_readvariableop_resource:
��>
/block_0_layer_0_biasadd_readvariableop_resource:	�B
.block_1_layer_0_matmul_readvariableop_resource:
��>
/block_1_layer_0_biasadd_readvariableop_resource:	�B
.block_2_layer_0_matmul_readvariableop_resource:
��>
/block_2_layer_0_biasadd_readvariableop_resource:	�B
.block_3_layer_0_matmul_readvariableop_resource:
��>
/block_3_layer_0_biasadd_readvariableop_resource:	�B
.block_4_layer_0_matmul_readvariableop_resource:
��>
/block_4_layer_0_biasadd_readvariableop_resource:	�B
.block_5_layer_0_matmul_readvariableop_resource:
��>
/block_5_layer_0_biasadd_readvariableop_resource:	�>
+dense_output_matmul_readvariableop_resource:	�:
,dense_output_biasadd_readvariableop_resource:
identity��&block_0_layer_0/BiasAdd/ReadVariableOp�%block_0_layer_0/MatMul/ReadVariableOp�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_1_layer_0/BiasAdd/ReadVariableOp�%block_1_layer_0/MatMul/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_2_layer_0/BiasAdd/ReadVariableOp�%block_2_layer_0/MatMul/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_3_layer_0/BiasAdd/ReadVariableOp�%block_3_layer_0/MatMul/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_4_layer_0/BiasAdd/ReadVariableOp�%block_4_layer_0/MatMul/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_5_layer_0/BiasAdd/ReadVariableOp�%block_5_layer_0/MatMul/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�)decorrelation_layer/MatMul/ReadVariableOp�#dense_output/BiasAdd/ReadVariableOp�"dense_output/MatMul/ReadVariableOp�5dense_output/kernel/Regularizer/Square/ReadVariableOp�"layer_input/BiasAdd/ReadVariableOp�!layer_input/MatMul/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�#mean_shift_layer/sub/ReadVariableOp�
#mean_shift_layer/sub/ReadVariableOpReadVariableOp,mean_shift_layer_sub_readvariableop_resource*
_output_shapes
:	*
dtype0�
mean_shift_layer/subSubinputs+mean_shift_layer/sub/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
)decorrelation_layer/MatMul/ReadVariableOpReadVariableOp2decorrelation_layer_matmul_readvariableop_resource*
_output_shapes

:		*
dtype0�
decorrelation_layer/MatMulMatMulmean_shift_layer/sub:z:01decorrelation_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
!layer_input/MatMul/ReadVariableOpReadVariableOp*layer_input_matmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
layer_input/MatMulMatMul$decorrelation_layer/MatMul:product:0)layer_input/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
"layer_input/BiasAdd/ReadVariableOpReadVariableOp+layer_input_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
layer_input/BiasAddBiasAddlayer_input/MatMul:product:0*layer_input/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
tf.math.softplus/SoftplusSoftpluslayer_input/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
%block_0_layer_0/MatMul/ReadVariableOpReadVariableOp.block_0_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_0_layer_0/MatMulMatMul'tf.math.softplus/Softplus:activations:0-block_0_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_0_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_0_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_0_layer_0/BiasAddBiasAdd block_0_layer_0/MatMul:product:0.block_0_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
add/addAddV2layer_input/BiasAdd:output:0 block_0_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������g
tf.math.softplus_1/SoftplusSoftplusadd/add:z:0*
T0*(
_output_shapes
:�����������
%block_1_layer_0/MatMul/ReadVariableOpReadVariableOp.block_1_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_1_layer_0/MatMulMatMul)tf.math.softplus_1/Softplus:activations:0-block_1_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_1_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_1_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_1_layer_0/BiasAddBiasAdd block_1_layer_0/MatMul:product:0.block_1_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������t
	add_1/addAddV2add/add:z:0 block_1_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������i
tf.math.softplus_2/SoftplusSoftplusadd_1/add:z:0*
T0*(
_output_shapes
:�����������
%block_2_layer_0/MatMul/ReadVariableOpReadVariableOp.block_2_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_2_layer_0/MatMulMatMul)tf.math.softplus_2/Softplus:activations:0-block_2_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_2_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_2_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_2_layer_0/BiasAddBiasAdd block_2_layer_0/MatMul:product:0.block_2_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
	add_2/addAddV2add_1/add:z:0 block_2_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������i
tf.math.softplus_3/SoftplusSoftplusadd_2/add:z:0*
T0*(
_output_shapes
:�����������
%block_3_layer_0/MatMul/ReadVariableOpReadVariableOp.block_3_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_3_layer_0/MatMulMatMul)tf.math.softplus_3/Softplus:activations:0-block_3_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_3_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_3_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_3_layer_0/BiasAddBiasAdd block_3_layer_0/MatMul:product:0.block_3_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
	add_3/addAddV2add_2/add:z:0 block_3_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������i
tf.math.softplus_4/SoftplusSoftplusadd_3/add:z:0*
T0*(
_output_shapes
:�����������
%block_4_layer_0/MatMul/ReadVariableOpReadVariableOp.block_4_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_4_layer_0/MatMulMatMul)tf.math.softplus_4/Softplus:activations:0-block_4_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_4_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_4_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_4_layer_0/BiasAddBiasAdd block_4_layer_0/MatMul:product:0.block_4_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
	add_4/addAddV2add_3/add:z:0 block_4_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������i
tf.math.softplus_5/SoftplusSoftplusadd_4/add:z:0*
T0*(
_output_shapes
:�����������
%block_5_layer_0/MatMul/ReadVariableOpReadVariableOp.block_5_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_5_layer_0/MatMulMatMul)tf.math.softplus_5/Softplus:activations:0-block_5_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_5_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_5_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_5_layer_0/BiasAddBiasAdd block_5_layer_0/MatMul:product:0.block_5_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
	add_5/addAddV2add_4/add:z:0 block_5_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
"dense_output/MatMul/ReadVariableOpReadVariableOp+dense_output_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_output/MatMulMatMuladd_5/add:z:0*dense_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
#dense_output/BiasAdd/ReadVariableOpReadVariableOp,dense_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_output/BiasAddBiasAdddense_output/MatMul:product:0+dense_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOp*layer_input_matmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOp+layer_input_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_0_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_0_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_1_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_1_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_2_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_2_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_3_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_3_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_4_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_4_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_5_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_5_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+dense_output_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: l
IdentityIdentitydense_output/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp'^block_0_layer_0/BiasAdd/ReadVariableOp&^block_0_layer_0/MatMul/ReadVariableOp7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_1_layer_0/BiasAdd/ReadVariableOp&^block_1_layer_0/MatMul/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_2_layer_0/BiasAdd/ReadVariableOp&^block_2_layer_0/MatMul/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_3_layer_0/BiasAdd/ReadVariableOp&^block_3_layer_0/MatMul/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_4_layer_0/BiasAdd/ReadVariableOp&^block_4_layer_0/MatMul/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_5_layer_0/BiasAdd/ReadVariableOp&^block_5_layer_0/MatMul/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp*^decorrelation_layer/MatMul/ReadVariableOp$^dense_output/BiasAdd/ReadVariableOp#^dense_output/MatMul/ReadVariableOp6^dense_output/kernel/Regularizer/Square/ReadVariableOp#^layer_input/BiasAdd/ReadVariableOp"^layer_input/MatMul/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp$^mean_shift_layer/sub/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 2P
&block_0_layer_0/BiasAdd/ReadVariableOp&block_0_layer_0/BiasAdd/ReadVariableOp2N
%block_0_layer_0/MatMul/ReadVariableOp%block_0_layer_0/MatMul/ReadVariableOp2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_1_layer_0/BiasAdd/ReadVariableOp&block_1_layer_0/BiasAdd/ReadVariableOp2N
%block_1_layer_0/MatMul/ReadVariableOp%block_1_layer_0/MatMul/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_2_layer_0/BiasAdd/ReadVariableOp&block_2_layer_0/BiasAdd/ReadVariableOp2N
%block_2_layer_0/MatMul/ReadVariableOp%block_2_layer_0/MatMul/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_3_layer_0/BiasAdd/ReadVariableOp&block_3_layer_0/BiasAdd/ReadVariableOp2N
%block_3_layer_0/MatMul/ReadVariableOp%block_3_layer_0/MatMul/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_4_layer_0/BiasAdd/ReadVariableOp&block_4_layer_0/BiasAdd/ReadVariableOp2N
%block_4_layer_0/MatMul/ReadVariableOp%block_4_layer_0/MatMul/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_5_layer_0/BiasAdd/ReadVariableOp&block_5_layer_0/BiasAdd/ReadVariableOp2N
%block_5_layer_0/MatMul/ReadVariableOp%block_5_layer_0/MatMul/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2V
)decorrelation_layer/MatMul/ReadVariableOp)decorrelation_layer/MatMul/ReadVariableOp2J
#dense_output/BiasAdd/ReadVariableOp#dense_output/BiasAdd/ReadVariableOp2H
"dense_output/MatMul/ReadVariableOp"dense_output/MatMul/ReadVariableOp2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2H
"layer_input/BiasAdd/ReadVariableOp"layer_input/BiasAdd/ReadVariableOp2F
!layer_input/MatMul/ReadVariableOp!layer_input/MatMul/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp2J
#mean_shift_layer/sub/ReadVariableOp#mean_shift_layer/sub/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
a__inference___backward_block_3_layer_0_layer_call_and_return_conditional_losses_93976592_93976605
placeholder6
2gradients_matmul_grad_matmul_matmul_readvariableop)
%gradients_matmul_grad_matmul_1_inputs
identity

identity_1

identity_2_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������u
"gradients/BiasAdd_grad/BiasAddGradBiasAddGradgradients/grad_ys_0:output:0*
T0*
_output_shapes	
:��
gradients/MatMul_grad/MatMulMatMulgradients/grad_ys_0:output:02gradients_matmul_grad_matmul_matmul_readvariableop*
T0*(
_output_shapes
:����������*
transpose_b(�
gradients/MatMul_grad/MatMul_1MatMul%gradients_matmul_grad_matmul_1_inputsgradients/grad_ys_0:output:0*
T0* 
_output_shapes
:
��*
transpose_a(o
IdentityIdentity&gradients/MatMul_grad/MatMul:product:0*
T0*(
_output_shapes
:����������k

Identity_1Identity(gradients/MatMul_grad/MatMul_1:product:0*
T0* 
_output_shapes
:
��i

Identity_2Identity+gradients/BiasAdd_grad/BiasAddGrad:output:0*
T0*
_output_shapes	
:�"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*G
_input_shapes6
4:����������:
��:����������*f
forward_function_nameMK__forward_block_3_layer_0_layer_call_and_return_conditional_losses_93976604:. *
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������
�
�
9__inference_ResNet_entropy_closure_layer_call_fn_93979630

inputs
unknown:	
	unknown_0:		
	unknown_1:		�
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:
��

unknown_12:	�

unknown_13:
��

unknown_14:	�

unknown_15:	�

unknown_16:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *]
fXRV
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976014o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
__inference_loss_fn_1_93980404J
;layer_input_bias_regularizer_square_readvariableop_resource:	�
identity��2layer_input/bias/Regularizer/Square/ReadVariableOp�
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOp;layer_input_bias_regularizer_square_readvariableop_resource*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: b
IdentityIdentity$layer_input/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: {
NoOpNoOp3^layer_input/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp
�
o
C__inference_add_2_layer_call_and_return_conditional_losses_93980186
inputs_0
inputs_1
identityS
addAddV2inputs_0inputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
�
�
A__forward_add_5_layer_call_and_return_conditional_losses_93976486
inputs_0

inputs_1_0
identity

inputs
inputs_1U
addAddV2inputs_0
inputs_1_0*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"
inputsinputs_0"
inputs_1
inputs_1_0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������*s
backward_function_nameYW__inference___backward_add_5_layer_call_and_return_conditional_losses_93976470_93976487:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_8_93980481U
Ablock_3_layer_0_kernel_regularizer_square_readvariableop_resource:
��
identity��8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpAblock_3_layer_0_kernel_regularizer_square_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: h
IdentityIdentity*block_3_layer_0/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp
�
�
A__forward_add_3_layer_call_and_return_conditional_losses_93976584
inputs_0

inputs_1_0
identity

inputs
inputs_1U
addAddV2inputs_0
inputs_1_0*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"
inputsinputs_0"
inputs_1
inputs_1_0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������*s
backward_function_nameYW__inference___backward_add_3_layer_call_and_return_conditional_losses_93976568_93976585:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_2_93980415U
Ablock_0_layer_0_kernel_regularizer_square_readvariableop_resource:
��
identity��8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpAblock_0_layer_0_kernel_regularizer_square_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: h
IdentityIdentity*block_0_layer_0/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp
�
�
__inference_loss_fn_12_93980525U
Ablock_5_layer_0_kernel_regularizer_square_readvariableop_resource:
��
identity��8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpAblock_5_layer_0_kernel_regularizer_square_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: h
IdentityIdentity*block_5_layer_0/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp
��
�-
$__inference__traced_restore_93981035
file_prefix$
assignvariableop_adam_iter:	 (
assignvariableop_1_adam_beta_1: (
assignvariableop_2_adam_beta_2: '
assignvariableop_3_adam_decay: /
%assignvariableop_4_adam_learning_rate: )
assignvariableop_5_variable:	/
assignvariableop_6_variable_1:		8
%assignvariableop_7_layer_input_kernel:		�2
#assignvariableop_8_layer_input_bias:	�=
)assignvariableop_9_block_0_layer_0_kernel:
��7
(assignvariableop_10_block_0_layer_0_bias:	�>
*assignvariableop_11_block_1_layer_0_kernel:
��7
(assignvariableop_12_block_1_layer_0_bias:	�>
*assignvariableop_13_block_2_layer_0_kernel:
��7
(assignvariableop_14_block_2_layer_0_bias:	�>
*assignvariableop_15_block_3_layer_0_kernel:
��7
(assignvariableop_16_block_3_layer_0_bias:	�>
*assignvariableop_17_block_4_layer_0_kernel:
��7
(assignvariableop_18_block_4_layer_0_bias:	�>
*assignvariableop_19_block_5_layer_0_kernel:
��7
(assignvariableop_20_block_5_layer_0_bias:	�:
'assignvariableop_21_dense_output_kernel:	�3
%assignvariableop_22_dense_output_bias:#
assignvariableop_23_total: #
assignvariableop_24_count: %
assignvariableop_25_total_1: %
assignvariableop_26_count_1: %
assignvariableop_27_total_2: %
assignvariableop_28_count_2: %
assignvariableop_29_total_3: %
assignvariableop_30_count_3: %
assignvariableop_31_total_4: %
assignvariableop_32_count_4: %
assignvariableop_33_total_5: %
assignvariableop_34_count_5: %
assignvariableop_35_total_6: %
assignvariableop_36_count_6: %
assignvariableop_37_total_7: %
assignvariableop_38_count_7: %
assignvariableop_39_total_8: %
assignvariableop_40_count_8: %
assignvariableop_41_total_9: %
assignvariableop_42_count_9: @
-assignvariableop_43_adam_layer_input_kernel_m:		�:
+assignvariableop_44_adam_layer_input_bias_m:	�E
1assignvariableop_45_adam_block_0_layer_0_kernel_m:
��>
/assignvariableop_46_adam_block_0_layer_0_bias_m:	�E
1assignvariableop_47_adam_block_1_layer_0_kernel_m:
��>
/assignvariableop_48_adam_block_1_layer_0_bias_m:	�E
1assignvariableop_49_adam_block_2_layer_0_kernel_m:
��>
/assignvariableop_50_adam_block_2_layer_0_bias_m:	�E
1assignvariableop_51_adam_block_3_layer_0_kernel_m:
��>
/assignvariableop_52_adam_block_3_layer_0_bias_m:	�E
1assignvariableop_53_adam_block_4_layer_0_kernel_m:
��>
/assignvariableop_54_adam_block_4_layer_0_bias_m:	�E
1assignvariableop_55_adam_block_5_layer_0_kernel_m:
��>
/assignvariableop_56_adam_block_5_layer_0_bias_m:	�A
.assignvariableop_57_adam_dense_output_kernel_m:	�:
,assignvariableop_58_adam_dense_output_bias_m:@
-assignvariableop_59_adam_layer_input_kernel_v:		�:
+assignvariableop_60_adam_layer_input_bias_v:	�E
1assignvariableop_61_adam_block_0_layer_0_kernel_v:
��>
/assignvariableop_62_adam_block_0_layer_0_bias_v:	�E
1assignvariableop_63_adam_block_1_layer_0_kernel_v:
��>
/assignvariableop_64_adam_block_1_layer_0_bias_v:	�E
1assignvariableop_65_adam_block_2_layer_0_kernel_v:
��>
/assignvariableop_66_adam_block_2_layer_0_bias_v:	�E
1assignvariableop_67_adam_block_3_layer_0_kernel_v:
��>
/assignvariableop_68_adam_block_3_layer_0_bias_v:	�E
1assignvariableop_69_adam_block_4_layer_0_kernel_v:
��>
/assignvariableop_70_adam_block_4_layer_0_bias_v:	�E
1assignvariableop_71_adam_block_5_layer_0_kernel_v:
��>
/assignvariableop_72_adam_block_5_layer_0_bias_v:	�A
.assignvariableop_73_adam_dense_output_kernel_v:	�:
,assignvariableop_74_adam_dense_output_bias_v:
identity_76��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_40�AssignVariableOp_41�AssignVariableOp_42�AssignVariableOp_43�AssignVariableOp_44�AssignVariableOp_45�AssignVariableOp_46�AssignVariableOp_47�AssignVariableOp_48�AssignVariableOp_49�AssignVariableOp_5�AssignVariableOp_50�AssignVariableOp_51�AssignVariableOp_52�AssignVariableOp_53�AssignVariableOp_54�AssignVariableOp_55�AssignVariableOp_56�AssignVariableOp_57�AssignVariableOp_58�AssignVariableOp_59�AssignVariableOp_6�AssignVariableOp_60�AssignVariableOp_61�AssignVariableOp_62�AssignVariableOp_63�AssignVariableOp_64�AssignVariableOp_65�AssignVariableOp_66�AssignVariableOp_67�AssignVariableOp_68�AssignVariableOp_69�AssignVariableOp_7�AssignVariableOp_70�AssignVariableOp_71�AssignVariableOp_72�AssignVariableOp_73�AssignVariableOp_74�AssignVariableOp_8�AssignVariableOp_9�"
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:L*
dtype0*�!
value�!B�!LB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB'variables/14/.ATTRIBUTES/VARIABLE_VALUEB'variables/15/.ATTRIBUTES/VARIABLE_VALUEB'variables/16/.ATTRIBUTES/VARIABLE_VALUEB'variables/17/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/3/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/3/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/4/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/4/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/5/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/5/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/6/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/6/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/7/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/7/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/8/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/8/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/9/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/9/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/16/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/17/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/16/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/17/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:L*
dtype0*�
value�B�LB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*Z
dtypesP
N2L	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOpAssignVariableOpassignvariableop_adam_iterIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOpassignvariableop_1_adam_beta_1Identity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOpassignvariableop_2_adam_beta_2Identity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOpassignvariableop_3_adam_decayIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOp%assignvariableop_4_adam_learning_rateIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpassignvariableop_5_variableIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOpassignvariableop_6_variable_1Identity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOp%assignvariableop_7_layer_input_kernelIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOp#assignvariableop_8_layer_input_biasIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOp)assignvariableop_9_block_0_layer_0_kernelIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOp(assignvariableop_10_block_0_layer_0_biasIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOp*assignvariableop_11_block_1_layer_0_kernelIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOp(assignvariableop_12_block_1_layer_0_biasIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOp*assignvariableop_13_block_2_layer_0_kernelIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_14AssignVariableOp(assignvariableop_14_block_2_layer_0_biasIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOp*assignvariableop_15_block_3_layer_0_kernelIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOp(assignvariableop_16_block_3_layer_0_biasIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOp*assignvariableop_17_block_4_layer_0_kernelIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOp(assignvariableop_18_block_4_layer_0_biasIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOp*assignvariableop_19_block_5_layer_0_kernelIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOp(assignvariableop_20_block_5_layer_0_biasIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOp'assignvariableop_21_dense_output_kernelIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOp%assignvariableop_22_dense_output_biasIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOpassignvariableop_23_totalIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOpassignvariableop_24_countIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOpassignvariableop_25_total_1Identity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOpassignvariableop_26_count_1Identity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOpassignvariableop_27_total_2Identity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOpassignvariableop_28_count_2Identity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_29AssignVariableOpassignvariableop_29_total_3Identity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_30AssignVariableOpassignvariableop_30_count_3Identity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_31AssignVariableOpassignvariableop_31_total_4Identity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_32AssignVariableOpassignvariableop_32_count_4Identity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_33AssignVariableOpassignvariableop_33_total_5Identity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_34AssignVariableOpassignvariableop_34_count_5Identity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_35AssignVariableOpassignvariableop_35_total_6Identity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_36AssignVariableOpassignvariableop_36_count_6Identity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_37AssignVariableOpassignvariableop_37_total_7Identity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_38AssignVariableOpassignvariableop_38_count_7Identity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_39AssignVariableOpassignvariableop_39_total_8Identity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_40AssignVariableOpassignvariableop_40_count_8Identity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_41AssignVariableOpassignvariableop_41_total_9Identity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_42AssignVariableOpassignvariableop_42_count_9Identity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_43AssignVariableOp-assignvariableop_43_adam_layer_input_kernel_mIdentity_43:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_44AssignVariableOp+assignvariableop_44_adam_layer_input_bias_mIdentity_44:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_45AssignVariableOp1assignvariableop_45_adam_block_0_layer_0_kernel_mIdentity_45:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_46AssignVariableOp/assignvariableop_46_adam_block_0_layer_0_bias_mIdentity_46:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_47AssignVariableOp1assignvariableop_47_adam_block_1_layer_0_kernel_mIdentity_47:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_48AssignVariableOp/assignvariableop_48_adam_block_1_layer_0_bias_mIdentity_48:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_49AssignVariableOp1assignvariableop_49_adam_block_2_layer_0_kernel_mIdentity_49:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_50IdentityRestoreV2:tensors:50"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_50AssignVariableOp/assignvariableop_50_adam_block_2_layer_0_bias_mIdentity_50:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_51IdentityRestoreV2:tensors:51"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_51AssignVariableOp1assignvariableop_51_adam_block_3_layer_0_kernel_mIdentity_51:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_52IdentityRestoreV2:tensors:52"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_52AssignVariableOp/assignvariableop_52_adam_block_3_layer_0_bias_mIdentity_52:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_53IdentityRestoreV2:tensors:53"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_53AssignVariableOp1assignvariableop_53_adam_block_4_layer_0_kernel_mIdentity_53:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_54IdentityRestoreV2:tensors:54"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_54AssignVariableOp/assignvariableop_54_adam_block_4_layer_0_bias_mIdentity_54:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_55IdentityRestoreV2:tensors:55"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_55AssignVariableOp1assignvariableop_55_adam_block_5_layer_0_kernel_mIdentity_55:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_56IdentityRestoreV2:tensors:56"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_56AssignVariableOp/assignvariableop_56_adam_block_5_layer_0_bias_mIdentity_56:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_57IdentityRestoreV2:tensors:57"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_57AssignVariableOp.assignvariableop_57_adam_dense_output_kernel_mIdentity_57:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_58IdentityRestoreV2:tensors:58"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_58AssignVariableOp,assignvariableop_58_adam_dense_output_bias_mIdentity_58:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_59IdentityRestoreV2:tensors:59"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_59AssignVariableOp-assignvariableop_59_adam_layer_input_kernel_vIdentity_59:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_60IdentityRestoreV2:tensors:60"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_60AssignVariableOp+assignvariableop_60_adam_layer_input_bias_vIdentity_60:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_61IdentityRestoreV2:tensors:61"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_61AssignVariableOp1assignvariableop_61_adam_block_0_layer_0_kernel_vIdentity_61:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_62IdentityRestoreV2:tensors:62"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_62AssignVariableOp/assignvariableop_62_adam_block_0_layer_0_bias_vIdentity_62:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_63IdentityRestoreV2:tensors:63"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_63AssignVariableOp1assignvariableop_63_adam_block_1_layer_0_kernel_vIdentity_63:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_64IdentityRestoreV2:tensors:64"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_64AssignVariableOp/assignvariableop_64_adam_block_1_layer_0_bias_vIdentity_64:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_65IdentityRestoreV2:tensors:65"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_65AssignVariableOp1assignvariableop_65_adam_block_2_layer_0_kernel_vIdentity_65:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_66IdentityRestoreV2:tensors:66"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_66AssignVariableOp/assignvariableop_66_adam_block_2_layer_0_bias_vIdentity_66:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_67IdentityRestoreV2:tensors:67"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_67AssignVariableOp1assignvariableop_67_adam_block_3_layer_0_kernel_vIdentity_67:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_68IdentityRestoreV2:tensors:68"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_68AssignVariableOp/assignvariableop_68_adam_block_3_layer_0_bias_vIdentity_68:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_69IdentityRestoreV2:tensors:69"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_69AssignVariableOp1assignvariableop_69_adam_block_4_layer_0_kernel_vIdentity_69:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_70IdentityRestoreV2:tensors:70"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_70AssignVariableOp/assignvariableop_70_adam_block_4_layer_0_bias_vIdentity_70:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_71IdentityRestoreV2:tensors:71"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_71AssignVariableOp1assignvariableop_71_adam_block_5_layer_0_kernel_vIdentity_71:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_72IdentityRestoreV2:tensors:72"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_72AssignVariableOp/assignvariableop_72_adam_block_5_layer_0_bias_vIdentity_72:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_73IdentityRestoreV2:tensors:73"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_73AssignVariableOp.assignvariableop_73_adam_dense_output_kernel_vIdentity_73:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_74IdentityRestoreV2:tensors:74"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_74AssignVariableOp,assignvariableop_74_adam_dense_output_bias_vIdentity_74:output:0"/device:CPU:0*
_output_shapes
 *
dtype01
NoOpNoOp"/device:CPU:0*
_output_shapes
 �
Identity_75Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_62^AssignVariableOp_63^AssignVariableOp_64^AssignVariableOp_65^AssignVariableOp_66^AssignVariableOp_67^AssignVariableOp_68^AssignVariableOp_69^AssignVariableOp_7^AssignVariableOp_70^AssignVariableOp_71^AssignVariableOp_72^AssignVariableOp_73^AssignVariableOp_74^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_76IdentityIdentity_75:output:0^NoOp_1*
T0*
_output_shapes
: �
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_62^AssignVariableOp_63^AssignVariableOp_64^AssignVariableOp_65^AssignVariableOp_66^AssignVariableOp_67^AssignVariableOp_68^AssignVariableOp_69^AssignVariableOp_7^AssignVariableOp_70^AssignVariableOp_71^AssignVariableOp_72^AssignVariableOp_73^AssignVariableOp_74^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 "#
identity_76Identity_76:output:0*�
_input_shapes�
�: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_5AssignVariableOp_52*
AssignVariableOp_50AssignVariableOp_502*
AssignVariableOp_51AssignVariableOp_512*
AssignVariableOp_52AssignVariableOp_522*
AssignVariableOp_53AssignVariableOp_532*
AssignVariableOp_54AssignVariableOp_542*
AssignVariableOp_55AssignVariableOp_552*
AssignVariableOp_56AssignVariableOp_562*
AssignVariableOp_57AssignVariableOp_572*
AssignVariableOp_58AssignVariableOp_582*
AssignVariableOp_59AssignVariableOp_592(
AssignVariableOp_6AssignVariableOp_62*
AssignVariableOp_60AssignVariableOp_602*
AssignVariableOp_61AssignVariableOp_612*
AssignVariableOp_62AssignVariableOp_622*
AssignVariableOp_63AssignVariableOp_632*
AssignVariableOp_64AssignVariableOp_642*
AssignVariableOp_65AssignVariableOp_652*
AssignVariableOp_66AssignVariableOp_662*
AssignVariableOp_67AssignVariableOp_672*
AssignVariableOp_68AssignVariableOp_682*
AssignVariableOp_69AssignVariableOp_692(
AssignVariableOp_7AssignVariableOp_72*
AssignVariableOp_70AssignVariableOp_702*
AssignVariableOp_71AssignVariableOp_712*
AssignVariableOp_72AssignVariableOp_722*
AssignVariableOp_73AssignVariableOp_732*
AssignVariableOp_74AssignVariableOp_742(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�
�
U__inference___backward_add_layer_call_and_return_conditional_losses_93976715_93976732
placeholder#
gradients_add_grad_shape_inputs'
#gradients_add_grad_shape_1_inputs_1
identity

identity_1_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������g
gradients/add_grad/ShapeShapegradients_add_grad_shape_inputs*
T0*
_output_shapes
:m
gradients/add_grad/Shape_1Shape#gradients_add_grad_shape_1_inputs_1*
T0*
_output_shapes
:�
(gradients/add_grad/BroadcastGradientArgsBroadcastGradientArgs!gradients/add_grad/Shape:output:0#gradients/add_grad/Shape_1:output:0*2
_output_shapes 
:���������:����������
gradients/add_grad/SumSumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
gradients/add_grad/ReshapeReshapegradients/add_grad/Sum:output:0!gradients/add_grad/Shape:output:0*
T0*(
_output_shapes
:�����������
gradients/add_grad/Sum_1Sumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
gradients/add_grad/Reshape_1Reshape!gradients/add_grad/Sum_1:output:0#gradients/add_grad/Shape_1:output:0*
T0*(
_output_shapes
:����������l
IdentityIdentity#gradients/add_grad/Reshape:output:0*
T0*(
_output_shapes
:����������p

Identity_1Identity%gradients/add_grad/Reshape_1:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������:����������:����������*Z
forward_function_nameA?__forward_add_layer_call_and_return_conditional_losses_93976731:. *
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������
�
�
M__inference_block_0_layer_0_layer_call_and_return_conditional_losses_93975330

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_6_93980459U
Ablock_2_layer_0_kernel_regularizer_square_readvariableop_resource:
��
identity��8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpAblock_2_layer_0_kernel_regularizer_square_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: h
IdentityIdentity*block_2_layer_0/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp
�
o
C__inference_add_4_layer_call_and_return_conditional_losses_93980296
inputs_0
inputs_1
identityS
addAddV2inputs_0inputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
�
�
__inference_loss_fn_5_93980448N
?block_1_layer_0_bias_regularizer_square_readvariableop_resource:	�
identity��6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp?block_1_layer_0_bias_regularizer_square_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: f
IdentityIdentity(block_1_layer_0/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 
NoOpNoOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp
�
�
6__inference_decorrelation_layer_layer_call_fn_93979971

inputs
unknown:		
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *Z
fURS
Q__inference_decorrelation_layer_layer_call_and_return_conditional_losses_93975275o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������	`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*(
_input_shapes
:���������	: 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
��
�
h__inference___backward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93977269_93977421
placeholder\
Xgradients_dense_output_statefulpartitionedcall_grad_dense_output_statefulpartitionedcall^
Zgradients_dense_output_statefulpartitionedcall_grad_dense_output_statefulpartitionedcall_1>
:gradients_add_5_partitionedcall_grad_add_5_partitionedcall@
<gradients_add_5_partitionedcall_grad_add_5_partitionedcall_1b
^gradients_block_5_layer_0_statefulpartitionedcall_grad_block_5_layer_0_statefulpartitionedcalld
`gradients_block_5_layer_0_statefulpartitionedcall_grad_block_5_layer_0_statefulpartitionedcall_1L
Hgradients_tf_math_softplus_5_softplus_grad_sigmoid_add_4_partitionedcall>
:gradients_add_4_partitionedcall_grad_add_4_partitionedcall@
<gradients_add_4_partitionedcall_grad_add_4_partitionedcall_1b
^gradients_block_4_layer_0_statefulpartitionedcall_grad_block_4_layer_0_statefulpartitionedcalld
`gradients_block_4_layer_0_statefulpartitionedcall_grad_block_4_layer_0_statefulpartitionedcall_1L
Hgradients_tf_math_softplus_4_softplus_grad_sigmoid_add_3_partitionedcall>
:gradients_add_3_partitionedcall_grad_add_3_partitionedcall@
<gradients_add_3_partitionedcall_grad_add_3_partitionedcall_1b
^gradients_block_3_layer_0_statefulpartitionedcall_grad_block_3_layer_0_statefulpartitionedcalld
`gradients_block_3_layer_0_statefulpartitionedcall_grad_block_3_layer_0_statefulpartitionedcall_1L
Hgradients_tf_math_softplus_3_softplus_grad_sigmoid_add_2_partitionedcall>
:gradients_add_2_partitionedcall_grad_add_2_partitionedcall@
<gradients_add_2_partitionedcall_grad_add_2_partitionedcall_1b
^gradients_block_2_layer_0_statefulpartitionedcall_grad_block_2_layer_0_statefulpartitionedcalld
`gradients_block_2_layer_0_statefulpartitionedcall_grad_block_2_layer_0_statefulpartitionedcall_1L
Hgradients_tf_math_softplus_2_softplus_grad_sigmoid_add_1_partitionedcall>
:gradients_add_1_partitionedcall_grad_add_1_partitionedcall@
<gradients_add_1_partitionedcall_grad_add_1_partitionedcall_1b
^gradients_block_1_layer_0_statefulpartitionedcall_grad_block_1_layer_0_statefulpartitionedcalld
`gradients_block_1_layer_0_statefulpartitionedcall_grad_block_1_layer_0_statefulpartitionedcall_1J
Fgradients_tf_math_softplus_1_softplus_grad_sigmoid_add_partitionedcall:
6gradients_add_partitionedcall_grad_add_partitionedcall<
8gradients_add_partitionedcall_grad_add_partitionedcall_1b
^gradients_block_0_layer_0_statefulpartitionedcall_grad_block_0_layer_0_statefulpartitionedcalld
`gradients_block_0_layer_0_statefulpartitionedcall_grad_block_0_layer_0_statefulpartitionedcall_1X
Tgradients_tf_math_softplus_softplus_grad_sigmoid_layer_input_statefulpartitionedcallZ
Vgradients_layer_input_statefulpartitionedcall_grad_layer_input_statefulpartitionedcall\
Xgradients_layer_input_statefulpartitionedcall_grad_layer_input_statefulpartitionedcall_1j
fgradients_decorrelation_layer_statefulpartitionedcall_grad_decorrelation_layer_statefulpartitionedcalll
hgradients_decorrelation_layer_statefulpartitionedcall_grad_decorrelation_layer_statefulpartitionedcall_1d
`gradients_mean_shift_layer_statefulpartitionedcall_grad_mean_shift_layer_statefulpartitionedcallf
bgradients_mean_shift_layer_statefulpartitionedcall_grad_mean_shift_layer_statefulpartitionedcall_1
identity

identity_1

identity_2

identity_3

identity_4

identity_5

identity_6

identity_7

identity_8

identity_9
identity_10
identity_11
identity_12
identity_13
identity_14
identity_15
identity_16
identity_17
identity_18^
gradients/grad_ys_0Identityplaceholder*
T0*'
_output_shapes
:����������
Cgradients/dense_output/StatefulPartitionedCall_grad/PartitionedCallPartitionedCallgradients/grad_ys_0:output:0Xgradients_dense_output_statefulpartitionedcall_grad_dense_output_statefulpartitionedcallZgradients_dense_output_statefulpartitionedcall_grad_dense_output_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *9
_output_shapes'
%:����������:	�:* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *g
fbR`
^__inference___backward_dense_output_layer_call_and_return_conditional_losses_93976449_93976462�
4gradients/add_5/PartitionedCall_grad/PartitionedCallPartitionedCallLgradients/dense_output/StatefulPartitionedCall_grad/PartitionedCall:output:0:gradients_add_5_partitionedcall_grad_add_5_partitionedcall<gradients_add_5_partitionedcall_grad_add_5_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_5_layer_call_and_return_conditional_losses_93976470_93976487�
Fgradients/block_5_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_5/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_5_layer_0_statefulpartitionedcall_grad_block_5_layer_0_statefulpartitionedcall`gradients_block_5_layer_0_statefulpartitionedcall_grad_block_5_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_5_layer_0_layer_call_and_return_conditional_losses_93976494_93976507�
2gradients/tf.math.softplus_5/Softplus_grad/SigmoidSigmoidHgradients_tf_math_softplus_5_softplus_grad_sigmoid_add_4_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_5/Softplus_grad/mulMulOgradients/block_5_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_5/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddNAddN=gradients/add_5/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_5/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_5/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
4gradients/add_4/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN:sum:0:gradients_add_4_partitionedcall_grad_add_4_partitionedcall<gradients_add_4_partitionedcall_grad_add_4_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_4_layer_call_and_return_conditional_losses_93976519_93976536�
Fgradients/block_4_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_4/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_4_layer_0_statefulpartitionedcall_grad_block_4_layer_0_statefulpartitionedcall`gradients_block_4_layer_0_statefulpartitionedcall_grad_block_4_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_4_layer_0_layer_call_and_return_conditional_losses_93976543_93976556�
2gradients/tf.math.softplus_4/Softplus_grad/SigmoidSigmoidHgradients_tf_math_softplus_4_softplus_grad_sigmoid_add_3_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_4/Softplus_grad/mulMulOgradients/block_4_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_4/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_1AddN=gradients/add_4/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_4/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_4/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
4gradients/add_3/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_1:sum:0:gradients_add_3_partitionedcall_grad_add_3_partitionedcall<gradients_add_3_partitionedcall_grad_add_3_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_3_layer_call_and_return_conditional_losses_93976568_93976585�
Fgradients/block_3_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_3/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_3_layer_0_statefulpartitionedcall_grad_block_3_layer_0_statefulpartitionedcall`gradients_block_3_layer_0_statefulpartitionedcall_grad_block_3_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_3_layer_0_layer_call_and_return_conditional_losses_93976592_93976605�
2gradients/tf.math.softplus_3/Softplus_grad/SigmoidSigmoidHgradients_tf_math_softplus_3_softplus_grad_sigmoid_add_2_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_3/Softplus_grad/mulMulOgradients/block_3_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_3/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_2AddN=gradients/add_3/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_3/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_3/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
4gradients/add_2/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_2:sum:0:gradients_add_2_partitionedcall_grad_add_2_partitionedcall<gradients_add_2_partitionedcall_grad_add_2_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_2_layer_call_and_return_conditional_losses_93976617_93976634�
Fgradients/block_2_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_2/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_2_layer_0_statefulpartitionedcall_grad_block_2_layer_0_statefulpartitionedcall`gradients_block_2_layer_0_statefulpartitionedcall_grad_block_2_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_2_layer_0_layer_call_and_return_conditional_losses_93976641_93976654�
2gradients/tf.math.softplus_2/Softplus_grad/SigmoidSigmoidHgradients_tf_math_softplus_2_softplus_grad_sigmoid_add_1_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_2/Softplus_grad/mulMulOgradients/block_2_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_2/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_3AddN=gradients/add_2/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_2/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_2/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
4gradients/add_1/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_3:sum:0:gradients_add_1_partitionedcall_grad_add_1_partitionedcall<gradients_add_1_partitionedcall_grad_add_1_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *`
f[RY
W__inference___backward_add_1_layer_call_and_return_conditional_losses_93976666_93976683�
Fgradients/block_1_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall=gradients/add_1/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_1_layer_0_statefulpartitionedcall_grad_block_1_layer_0_statefulpartitionedcall`gradients_block_1_layer_0_statefulpartitionedcall_grad_block_1_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_1_layer_0_layer_call_and_return_conditional_losses_93976690_93976703�
2gradients/tf.math.softplus_1/Softplus_grad/SigmoidSigmoidFgradients_tf_math_softplus_1_softplus_grad_sigmoid_add_partitionedcall*
T0*(
_output_shapes
:�����������
.gradients/tf.math.softplus_1/Softplus_grad/mulMulOgradients/block_1_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:06gradients/tf.math.softplus_1/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_4AddN=gradients/add_1/PartitionedCall_grad/PartitionedCall:output:02gradients/tf.math.softplus_1/Softplus_grad/mul:z:0*
N*
T0*G
_class=
;9loc:@gradients/add_1/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
2gradients/add/PartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_4:sum:06gradients_add_partitionedcall_grad_add_partitionedcall8gradients_add_partitionedcall_grad_add_partitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *<
_output_shapes*
(:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *^
fYRW
U__inference___backward_add_layer_call_and_return_conditional_losses_93976715_93976732�
Fgradients/block_0_layer_0/StatefulPartitionedCall_grad/PartitionedCallPartitionedCall;gradients/add/PartitionedCall_grad/PartitionedCall:output:1^gradients_block_0_layer_0_statefulpartitionedcall_grad_block_0_layer_0_statefulpartitionedcall`gradients_block_0_layer_0_statefulpartitionedcall_grad_block_0_layer_0_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *;
_output_shapes)
':����������:
��:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *j
feRc
a__inference___backward_block_0_layer_0_layer_call_and_return_conditional_losses_93976739_93976752�
0gradients/tf.math.softplus/Softplus_grad/SigmoidSigmoidTgradients_tf_math_softplus_softplus_grad_sigmoid_layer_input_statefulpartitionedcall*
T0*(
_output_shapes
:�����������
,gradients/tf.math.softplus/Softplus_grad/mulMulOgradients/block_0_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:04gradients/tf.math.softplus/Softplus_grad/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
gradients/AddN_5AddN;gradients/add/PartitionedCall_grad/PartitionedCall:output:00gradients/tf.math.softplus/Softplus_grad/mul:z:0*
N*
T0*E
_class;
97loc:@gradients/add/PartitionedCall_grad/PartitionedCall*(
_output_shapes
:�����������
Bgradients/layer_input/StatefulPartitionedCall_grad/PartitionedCallPartitionedCallgradients/AddN_5:sum:0Vgradients_layer_input_statefulpartitionedcall_grad_layer_input_statefulpartitionedcallXgradients_layer_input_statefulpartitionedcall_grad_layer_input_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *9
_output_shapes'
%:���������	:		�:�* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *f
faR_
]__inference___backward_layer_input_layer_call_and_return_conditional_losses_93976764_93976777�
Jgradients/decorrelation_layer/StatefulPartitionedCall_grad/PartitionedCallPartitionedCallKgradients/layer_input/StatefulPartitionedCall_grad/PartitionedCall:output:0fgradients_decorrelation_layer_statefulpartitionedcall_grad_decorrelation_layer_statefulpartitionedcallhgradients_decorrelation_layer_statefulpartitionedcall_grad_decorrelation_layer_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *1
_output_shapes
:���������	:		* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *n
fiRg
e__inference___backward_decorrelation_layer_layer_call_and_return_conditional_losses_93976785_93976796�
Ggradients/mean_shift_layer/StatefulPartitionedCall_grad/PartitionedCallPartitionedCallSgradients/decorrelation_layer/StatefulPartitionedCall_grad/PartitionedCall:output:0`gradients_mean_shift_layer_statefulpartitionedcall_grad_mean_shift_layer_statefulpartitionedcallbgradients_mean_shift_layer_statefulpartitionedcall_grad_mean_shift_layer_statefulpartitionedcall_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *-
_output_shapes
:���������	:	* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *k
ffRd
b__inference___backward_mean_shift_layer_layer_call_and_return_conditional_losses_93976803_93976821�
IdentityIdentityPgradients/mean_shift_layer/StatefulPartitionedCall_grad/PartitionedCall:output:0*
T0*'
_output_shapes
:���������	�

Identity_1IdentityPgradients/mean_shift_layer/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0*
_output_shapes
:	�

Identity_2IdentitySgradients/decorrelation_layer/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0*
_output_shapes

:		�

Identity_3IdentityKgradients/layer_input/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0*
_output_shapes
:		��

Identity_4IdentityKgradients/layer_input/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��

Identity_5IdentityOgradients/block_0_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���

Identity_6IdentityOgradients/block_0_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��

Identity_7IdentityOgradients/block_1_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���

Identity_8IdentityOgradients/block_1_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��

Identity_9IdentityOgradients/block_2_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���
Identity_10IdentityOgradients/block_2_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��
Identity_11IdentityOgradients/block_3_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���
Identity_12IdentityOgradients/block_3_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��
Identity_13IdentityOgradients/block_4_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���
Identity_14IdentityOgradients/block_4_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��
Identity_15IdentityOgradients/block_5_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0* 
_output_shapes
:
���
Identity_16IdentityOgradients/block_5_layer_0/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes	
:��
Identity_17IdentityLgradients/dense_output/StatefulPartitionedCall_grad/PartitionedCall:output:1*
T0*
_output_shapes
:	��
Identity_18IdentityLgradients/dense_output/StatefulPartitionedCall_grad/PartitionedCall:output:2*
T0*
_output_shapes
:"
identityIdentity:output:0"!

identity_1Identity_1:output:0"#
identity_10Identity_10:output:0"#
identity_11Identity_11:output:0"#
identity_12Identity_12:output:0"#
identity_13Identity_13:output:0"#
identity_14Identity_14:output:0"#
identity_15Identity_15:output:0"#
identity_16Identity_16:output:0"#
identity_17Identity_17:output:0"#
identity_18Identity_18:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0"!

identity_4Identity_4:output:0"!

identity_5Identity_5:output:0"!

identity_6Identity_6:output:0"!

identity_7Identity_7:output:0"!

identity_8Identity_8:output:0"!

identity_9Identity_9:output:0*(
_construction_contextkEagerRuntime*�
_input_shapes�
�:���������:	�:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:		�:���������	:		:���������	:���������	:	*m
forward_function_nameTR__forward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93977420:- )
'
_output_shapes
:���������:%!

_output_shapes
:	�:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.	*
(
_output_shapes
:����������:&
"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������:. *
(
_output_shapes
:����������:%!!

_output_shapes
:		�:-")
'
_output_shapes
:���������	:$# 

_output_shapes

:		:-$)
'
_output_shapes
:���������	:-%)
'
_output_shapes
:���������	: &

_output_shapes
:	
�
�
2__inference_block_2_layer_0_layer_call_fn_93980152

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_2_layer_0_layer_call_and_return_conditional_losses_93975404p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
m
C__inference_add_3_layer_call_and_return_conditional_losses_93975453

inputs
inputs_1
identityQ
addAddV2inputsinputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
9__inference_ResNet_entropy_closure_layer_call_fn_93976094
input_1
unknown:	
	unknown_0:		
	unknown_1:		�
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:
��

unknown_12:	�

unknown_13:
��

unknown_14:	�

unknown_15:	�

unknown_16:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *]
fXRV
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976014o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1
�
o
C__inference_add_3_layer_call_and_return_conditional_losses_93980241
inputs_0
inputs_1
identityS
addAddV2inputs_0inputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
�
�
0__inference_sobolev_model_layer_call_fn_93977173
input_1
unknown:	
	unknown_0:		
	unknown_1:		�
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:
��

unknown_12:	�

unknown_13:
��

unknown_14:	�

unknown_15:	�

unknown_16:

unknown_17

unknown_18

unknown_19
identity

identity_1

identity_2��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19*!
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������	:���������	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93977124o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������	q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������	`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
�
M__inference_block_3_layer_0_layer_call_and_return_conditional_losses_93980229

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93979790

inputs:
,mean_shift_layer_sub_readvariableop_resource:	D
2decorrelation_layer_matmul_readvariableop_resource:		=
*layer_input_matmul_readvariableop_resource:		�:
+layer_input_biasadd_readvariableop_resource:	�B
.block_0_layer_0_matmul_readvariableop_resource:
��>
/block_0_layer_0_biasadd_readvariableop_resource:	�B
.block_1_layer_0_matmul_readvariableop_resource:
��>
/block_1_layer_0_biasadd_readvariableop_resource:	�B
.block_2_layer_0_matmul_readvariableop_resource:
��>
/block_2_layer_0_biasadd_readvariableop_resource:	�B
.block_3_layer_0_matmul_readvariableop_resource:
��>
/block_3_layer_0_biasadd_readvariableop_resource:	�B
.block_4_layer_0_matmul_readvariableop_resource:
��>
/block_4_layer_0_biasadd_readvariableop_resource:	�B
.block_5_layer_0_matmul_readvariableop_resource:
��>
/block_5_layer_0_biasadd_readvariableop_resource:	�>
+dense_output_matmul_readvariableop_resource:	�:
,dense_output_biasadd_readvariableop_resource:
identity��&block_0_layer_0/BiasAdd/ReadVariableOp�%block_0_layer_0/MatMul/ReadVariableOp�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_1_layer_0/BiasAdd/ReadVariableOp�%block_1_layer_0/MatMul/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_2_layer_0/BiasAdd/ReadVariableOp�%block_2_layer_0/MatMul/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_3_layer_0/BiasAdd/ReadVariableOp�%block_3_layer_0/MatMul/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_4_layer_0/BiasAdd/ReadVariableOp�%block_4_layer_0/MatMul/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�&block_5_layer_0/BiasAdd/ReadVariableOp�%block_5_layer_0/MatMul/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�)decorrelation_layer/MatMul/ReadVariableOp�#dense_output/BiasAdd/ReadVariableOp�"dense_output/MatMul/ReadVariableOp�5dense_output/kernel/Regularizer/Square/ReadVariableOp�"layer_input/BiasAdd/ReadVariableOp�!layer_input/MatMul/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�#mean_shift_layer/sub/ReadVariableOp�
#mean_shift_layer/sub/ReadVariableOpReadVariableOp,mean_shift_layer_sub_readvariableop_resource*
_output_shapes
:	*
dtype0�
mean_shift_layer/subSubinputs+mean_shift_layer/sub/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
)decorrelation_layer/MatMul/ReadVariableOpReadVariableOp2decorrelation_layer_matmul_readvariableop_resource*
_output_shapes

:		*
dtype0�
decorrelation_layer/MatMulMatMulmean_shift_layer/sub:z:01decorrelation_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
!layer_input/MatMul/ReadVariableOpReadVariableOp*layer_input_matmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
layer_input/MatMulMatMul$decorrelation_layer/MatMul:product:0)layer_input/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
"layer_input/BiasAdd/ReadVariableOpReadVariableOp+layer_input_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
layer_input/BiasAddBiasAddlayer_input/MatMul:product:0*layer_input/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
tf.math.softplus/SoftplusSoftpluslayer_input/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
%block_0_layer_0/MatMul/ReadVariableOpReadVariableOp.block_0_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_0_layer_0/MatMulMatMul'tf.math.softplus/Softplus:activations:0-block_0_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_0_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_0_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_0_layer_0/BiasAddBiasAdd block_0_layer_0/MatMul:product:0.block_0_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
add/addAddV2layer_input/BiasAdd:output:0 block_0_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������g
tf.math.softplus_1/SoftplusSoftplusadd/add:z:0*
T0*(
_output_shapes
:�����������
%block_1_layer_0/MatMul/ReadVariableOpReadVariableOp.block_1_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_1_layer_0/MatMulMatMul)tf.math.softplus_1/Softplus:activations:0-block_1_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_1_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_1_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_1_layer_0/BiasAddBiasAdd block_1_layer_0/MatMul:product:0.block_1_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������t
	add_1/addAddV2add/add:z:0 block_1_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������i
tf.math.softplus_2/SoftplusSoftplusadd_1/add:z:0*
T0*(
_output_shapes
:�����������
%block_2_layer_0/MatMul/ReadVariableOpReadVariableOp.block_2_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_2_layer_0/MatMulMatMul)tf.math.softplus_2/Softplus:activations:0-block_2_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_2_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_2_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_2_layer_0/BiasAddBiasAdd block_2_layer_0/MatMul:product:0.block_2_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
	add_2/addAddV2add_1/add:z:0 block_2_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������i
tf.math.softplus_3/SoftplusSoftplusadd_2/add:z:0*
T0*(
_output_shapes
:�����������
%block_3_layer_0/MatMul/ReadVariableOpReadVariableOp.block_3_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_3_layer_0/MatMulMatMul)tf.math.softplus_3/Softplus:activations:0-block_3_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_3_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_3_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_3_layer_0/BiasAddBiasAdd block_3_layer_0/MatMul:product:0.block_3_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
	add_3/addAddV2add_2/add:z:0 block_3_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������i
tf.math.softplus_4/SoftplusSoftplusadd_3/add:z:0*
T0*(
_output_shapes
:�����������
%block_4_layer_0/MatMul/ReadVariableOpReadVariableOp.block_4_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_4_layer_0/MatMulMatMul)tf.math.softplus_4/Softplus:activations:0-block_4_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_4_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_4_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_4_layer_0/BiasAddBiasAdd block_4_layer_0/MatMul:product:0.block_4_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
	add_4/addAddV2add_3/add:z:0 block_4_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:����������i
tf.math.softplus_5/SoftplusSoftplusadd_4/add:z:0*
T0*(
_output_shapes
:�����������
%block_5_layer_0/MatMul/ReadVariableOpReadVariableOp.block_5_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
block_5_layer_0/MatMulMatMul)tf.math.softplus_5/Softplus:activations:0-block_5_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&block_5_layer_0/BiasAdd/ReadVariableOpReadVariableOp/block_5_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
block_5_layer_0/BiasAddBiasAdd block_5_layer_0/MatMul:product:0.block_5_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������v
	add_5/addAddV2add_4/add:z:0 block_5_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
"dense_output/MatMul/ReadVariableOpReadVariableOp+dense_output_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_output/MatMulMatMuladd_5/add:z:0*dense_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
#dense_output/BiasAdd/ReadVariableOpReadVariableOp,dense_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_output/BiasAddBiasAdddense_output/MatMul:product:0+dense_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOp*layer_input_matmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOp+layer_input_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_0_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_0_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_1_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_1_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_2_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_2_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_3_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_3_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_4_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_4_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOp.block_5_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp/block_5_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+dense_output_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: l
IdentityIdentitydense_output/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp'^block_0_layer_0/BiasAdd/ReadVariableOp&^block_0_layer_0/MatMul/ReadVariableOp7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_1_layer_0/BiasAdd/ReadVariableOp&^block_1_layer_0/MatMul/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_2_layer_0/BiasAdd/ReadVariableOp&^block_2_layer_0/MatMul/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_3_layer_0/BiasAdd/ReadVariableOp&^block_3_layer_0/MatMul/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_4_layer_0/BiasAdd/ReadVariableOp&^block_4_layer_0/MatMul/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp'^block_5_layer_0/BiasAdd/ReadVariableOp&^block_5_layer_0/MatMul/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp*^decorrelation_layer/MatMul/ReadVariableOp$^dense_output/BiasAdd/ReadVariableOp#^dense_output/MatMul/ReadVariableOp6^dense_output/kernel/Regularizer/Square/ReadVariableOp#^layer_input/BiasAdd/ReadVariableOp"^layer_input/MatMul/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp$^mean_shift_layer/sub/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 2P
&block_0_layer_0/BiasAdd/ReadVariableOp&block_0_layer_0/BiasAdd/ReadVariableOp2N
%block_0_layer_0/MatMul/ReadVariableOp%block_0_layer_0/MatMul/ReadVariableOp2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_1_layer_0/BiasAdd/ReadVariableOp&block_1_layer_0/BiasAdd/ReadVariableOp2N
%block_1_layer_0/MatMul/ReadVariableOp%block_1_layer_0/MatMul/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_2_layer_0/BiasAdd/ReadVariableOp&block_2_layer_0/BiasAdd/ReadVariableOp2N
%block_2_layer_0/MatMul/ReadVariableOp%block_2_layer_0/MatMul/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_3_layer_0/BiasAdd/ReadVariableOp&block_3_layer_0/BiasAdd/ReadVariableOp2N
%block_3_layer_0/MatMul/ReadVariableOp%block_3_layer_0/MatMul/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_4_layer_0/BiasAdd/ReadVariableOp&block_4_layer_0/BiasAdd/ReadVariableOp2N
%block_4_layer_0/MatMul/ReadVariableOp%block_4_layer_0/MatMul/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2P
&block_5_layer_0/BiasAdd/ReadVariableOp&block_5_layer_0/BiasAdd/ReadVariableOp2N
%block_5_layer_0/MatMul/ReadVariableOp%block_5_layer_0/MatMul/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2V
)decorrelation_layer/MatMul/ReadVariableOp)decorrelation_layer/MatMul/ReadVariableOp2J
#dense_output/BiasAdd/ReadVariableOp#dense_output/BiasAdd/ReadVariableOp2H
"dense_output/MatMul/ReadVariableOp"dense_output/MatMul/ReadVariableOp2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2H
"layer_input/BiasAdd/ReadVariableOp"layer_input/BiasAdd/ReadVariableOp2F
!layer_input/MatMul/ReadVariableOp!layer_input/MatMul/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp2J
#mean_shift_layer/sub/ReadVariableOp#mean_shift_layer/sub/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
2__inference_block_0_layer_0_layer_call_fn_93980042

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_0_layer_0_layer_call_and_return_conditional_losses_93975330p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
M__inference_block_1_layer_0_layer_call_and_return_conditional_losses_93975367

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_7_93980470N
?block_2_layer_0_bias_regularizer_square_readvariableop_resource:	�
identity��6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp?block_2_layer_0_bias_regularizer_square_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: f
IdentityIdentity(block_2_layer_0/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 
NoOpNoOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp
��
�
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93979458
xQ
Cresnet_entropy_closure_mean_shift_layer_sub_readvariableop_resource:	[
Iresnet_entropy_closure_decorrelation_layer_matmul_readvariableop_resource:		T
Aresnet_entropy_closure_layer_input_matmul_readvariableop_resource:		�Q
Bresnet_entropy_closure_layer_input_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_0_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_0_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_1_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_1_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_2_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_2_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_3_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_3_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_4_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_4_layer_0_biasadd_readvariableop_resource:	�Y
Eresnet_entropy_closure_block_5_layer_0_matmul_readvariableop_resource:
��U
Fresnet_entropy_closure_block_5_layer_0_biasadd_readvariableop_resource:	�U
Bresnet_entropy_closure_dense_output_matmul_readvariableop_resource:	�Q
Cresnet_entropy_closure_dense_output_biasadd_readvariableop_resource:
unknown
tensordot_1_b
mul_1_x
identity

identity_1

identity_2��=ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp�=ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp�<ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp�@ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp�:ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp�9ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp�9ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp�8ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp�:ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�checked�	checked_1�5dense_output/kernel/Regularizer/Square/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�
:ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOpReadVariableOpCresnet_entropy_closure_mean_shift_layer_sub_readvariableop_resource*
_output_shapes
:	*
dtype0�
+ResNet_entropy_closure/mean_shift_layer/subSubxBResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
@ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOpReadVariableOpIresnet_entropy_closure_decorrelation_layer_matmul_readvariableop_resource*
_output_shapes

:		*
dtype0�
1ResNet_entropy_closure/decorrelation_layer/MatMulMatMul/ResNet_entropy_closure/mean_shift_layer/sub:z:0HResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	�
8ResNet_entropy_closure/layer_input/MatMul/ReadVariableOpReadVariableOpAresnet_entropy_closure_layer_input_matmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
)ResNet_entropy_closure/layer_input/MatMulMatMul;ResNet_entropy_closure/decorrelation_layer/MatMul:product:0@ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
9ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOpReadVariableOpBresnet_entropy_closure_layer_input_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
*ResNet_entropy_closure/layer_input/BiasAddBiasAdd3ResNet_entropy_closure/layer_input/MatMul:product:0AResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
0ResNet_entropy_closure/tf.math.softplus/SoftplusSoftplus3ResNet_entropy_closure/layer_input/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_0_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_0_layer_0/MatMulMatMul>ResNet_entropy_closure/tf.math.softplus/Softplus:activations:0DResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_0_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_0_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_0_layer_0/MatMul:product:0EResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
ResNet_entropy_closure/add/addAddV23ResNet_entropy_closure/layer_input/BiasAdd:output:07ResNet_entropy_closure/block_0_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_1/SoftplusSoftplus"ResNet_entropy_closure/add/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_1_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_1_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_1/Softplus:activations:0DResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_1_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_1_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_1_layer_0/MatMul:product:0EResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_1/addAddV2"ResNet_entropy_closure/add/add:z:07ResNet_entropy_closure/block_1_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_2/SoftplusSoftplus$ResNet_entropy_closure/add_1/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_2_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_2_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_2/Softplus:activations:0DResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_2_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_2_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_2_layer_0/MatMul:product:0EResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_2/addAddV2$ResNet_entropy_closure/add_1/add:z:07ResNet_entropy_closure/block_2_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_3/SoftplusSoftplus$ResNet_entropy_closure/add_2/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_3_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_3_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_3/Softplus:activations:0DResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_3_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_3_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_3_layer_0/MatMul:product:0EResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_3/addAddV2$ResNet_entropy_closure/add_2/add:z:07ResNet_entropy_closure/block_3_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_4/SoftplusSoftplus$ResNet_entropy_closure/add_3/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_4_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_4_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_4/Softplus:activations:0DResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_4_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_4_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_4_layer_0/MatMul:product:0EResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_4/addAddV2$ResNet_entropy_closure/add_3/add:z:07ResNet_entropy_closure/block_4_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
2ResNet_entropy_closure/tf.math.softplus_5/SoftplusSoftplus$ResNet_entropy_closure/add_4/add:z:0*
T0*(
_output_shapes
:�����������
<ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_5_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
-ResNet_entropy_closure/block_5_layer_0/MatMulMatMul@ResNet_entropy_closure/tf.math.softplus_5/Softplus:activations:0DResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
=ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_5_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
.ResNet_entropy_closure/block_5_layer_0/BiasAddBiasAdd7ResNet_entropy_closure/block_5_layer_0/MatMul:product:0EResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
 ResNet_entropy_closure/add_5/addAddV2$ResNet_entropy_closure/add_4/add:z:07ResNet_entropy_closure/block_5_layer_0/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
9ResNet_entropy_closure/dense_output/MatMul/ReadVariableOpReadVariableOpBresnet_entropy_closure_dense_output_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
*ResNet_entropy_closure/dense_output/MatMulMatMul$ResNet_entropy_closure/add_5/add:z:0AResNet_entropy_closure/dense_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
:ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOpReadVariableOpCresnet_entropy_closure_dense_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
+ResNet_entropy_closure/dense_output/BiasAddBiasAdd4ResNet_entropy_closure/dense_output/MatMul:product:0BResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������s
ones_like/ShapeShape4ResNet_entropy_closure/dense_output/BiasAdd:output:0*
T0*
_output_shapes
:T
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?w
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:����������
Egradient_tape/ResNet_entropy_closure/dense_output/BiasAdd/BiasAddGradBiasAddGradones_like:output:0*
T0*
_output_shapes
:�
?gradient_tape/ResNet_entropy_closure/dense_output/MatMul/MatMulMatMulones_like:output:0AResNet_entropy_closure/dense_output/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Agradient_tape/ResNet_entropy_closure/dense_output/MatMul/MatMul_1MatMul$ResNet_entropy_closure/add_5/add:z:0ones_like:output:0*
T0*
_output_shapes
:	�*
transpose_a(�
4gradient_tape/ResNet_entropy_closure/add_5/add/ShapeShape$ResNet_entropy_closure/add_4/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_5/add/Shape_1Shape7ResNet_entropy_closure/block_5_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_5/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_5/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_5/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_5/add/SumSumIgradient_tape/ResNet_entropy_closure/dense_output/MatMul/MatMul:product:0Igradient_tape/ResNet_entropy_closure/add_5/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_5/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_5/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_5/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_5/add/Sum_1SumIgradient_tape/ResNet_entropy_closure/dense_output/MatMul/MatMul:product:0Igradient_tape/ResNet_entropy_closure/add_5/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_5/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_5/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_5/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_5_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_5/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_5_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_5/add/Reshape_1:output:0DResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_5_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_5/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_5/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_5/SigmoidSigmoid$ResNet_entropy_closure/add_4/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_5/mulMulLgradient_tape/ResNet_entropy_closure/block_5_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_5/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddNAddN?gradient_tape/ResNet_entropy_closure/add_5/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_5/mul:z:0*
N*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_4/add/ShapeShape$ResNet_entropy_closure/add_3/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_4/add/Shape_1Shape7ResNet_entropy_closure/block_4_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_4/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_4/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_4/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_4/add/SumSum
AddN:sum:0Igradient_tape/ResNet_entropy_closure/add_4/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_4/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_4/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_4/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_4/add/Sum_1Sum
AddN:sum:0Igradient_tape/ResNet_entropy_closure/add_4/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_4/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_4/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_4/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_4_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_4/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_4_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_4/add/Reshape_1:output:0DResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_4_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_4/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_4/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_4/SigmoidSigmoid$ResNet_entropy_closure/add_3/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_4/mulMulLgradient_tape/ResNet_entropy_closure/block_4_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_4/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_1AddN?gradient_tape/ResNet_entropy_closure/add_4/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_4/mul:z:0*
N*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_3/add/ShapeShape$ResNet_entropy_closure/add_2/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_3/add/Shape_1Shape7ResNet_entropy_closure/block_3_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_3/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_3/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_3/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_3/add/SumSumAddN_1:sum:0Igradient_tape/ResNet_entropy_closure/add_3/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_3/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_3/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_3/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_3/add/Sum_1SumAddN_1:sum:0Igradient_tape/ResNet_entropy_closure/add_3/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_3/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_3/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_3/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_3_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_3/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_3_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_3/add/Reshape_1:output:0DResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_3_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_3/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_3/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_3/SigmoidSigmoid$ResNet_entropy_closure/add_2/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_3/mulMulLgradient_tape/ResNet_entropy_closure/block_3_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_3/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_2AddN?gradient_tape/ResNet_entropy_closure/add_3/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_3/mul:z:0*
N*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_2/add/ShapeShape$ResNet_entropy_closure/add_1/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_2/add/Shape_1Shape7ResNet_entropy_closure/block_2_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_2/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_2/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_2/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_2/add/SumSumAddN_2:sum:0Igradient_tape/ResNet_entropy_closure/add_2/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_2/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_2/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_2/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_2/add/Sum_1SumAddN_2:sum:0Igradient_tape/ResNet_entropy_closure/add_2/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_2/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_2/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_2/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_2_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_2/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_2_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_2/add/Reshape_1:output:0DResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_2_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_2/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_2/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_2/SigmoidSigmoid$ResNet_entropy_closure/add_1/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_2/mulMulLgradient_tape/ResNet_entropy_closure/block_2_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_2/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_3AddN?gradient_tape/ResNet_entropy_closure/add_2/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_2/mul:z:0*
N*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_1/add/ShapeShape"ResNet_entropy_closure/add/add:z:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_1/add/Shape_1Shape7ResNet_entropy_closure/block_1_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Dgradient_tape/ResNet_entropy_closure/add_1/add/BroadcastGradientArgsBroadcastGradientArgs=gradient_tape/ResNet_entropy_closure/add_1/add/Shape:output:0?gradient_tape/ResNet_entropy_closure/add_1/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
2gradient_tape/ResNet_entropy_closure/add_1/add/SumSumAddN_3:sum:0Igradient_tape/ResNet_entropy_closure/add_1/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add_1/add/ReshapeReshape;gradient_tape/ResNet_entropy_closure/add_1/add/Sum:output:0=gradient_tape/ResNet_entropy_closure/add_1/add/Shape:output:0*
T0*(
_output_shapes
:�����������
4gradient_tape/ResNet_entropy_closure/add_1/add/Sum_1SumAddN_3:sum:0Igradient_tape/ResNet_entropy_closure/add_1/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
8gradient_tape/ResNet_entropy_closure/add_1/add/Reshape_1Reshape=gradient_tape/ResNet_entropy_closure/add_1/add/Sum_1:output:0?gradient_tape/ResNet_entropy_closure/add_1/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_1_layer_0/BiasAdd/BiasAddGradBiasAddGradAgradient_tape/ResNet_entropy_closure/add_1/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_1_layer_0/MatMul/MatMulMatMulAgradient_tape/ResNet_entropy_closure/add_1/add/Reshape_1:output:0DResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_1_layer_0/MatMul/MatMul_1MatMul@ResNet_entropy_closure/tf.math.softplus_1/Softplus:activations:0Agradient_tape/ResNet_entropy_closure/add_1/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
?gradient_tape/ResNet_entropy_closure/tf.math.softplus_1/SigmoidSigmoid"ResNet_entropy_closure/add/add:z:0*
T0*(
_output_shapes
:�����������
;gradient_tape/ResNet_entropy_closure/tf.math.softplus_1/mulMulLgradient_tape/ResNet_entropy_closure/block_1_layer_0/MatMul/MatMul:product:0Cgradient_tape/ResNet_entropy_closure/tf.math.softplus_1/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_4AddN?gradient_tape/ResNet_entropy_closure/add_1/add/Reshape:output:0?gradient_tape/ResNet_entropy_closure/tf.math.softplus_1/mul:z:0*
N*
T0*(
_output_shapes
:�����������
2gradient_tape/ResNet_entropy_closure/add/add/ShapeShape3ResNet_entropy_closure/layer_input/BiasAdd:output:0*
T0*
_output_shapes
:�
4gradient_tape/ResNet_entropy_closure/add/add/Shape_1Shape7ResNet_entropy_closure/block_0_layer_0/BiasAdd:output:0*
T0*
_output_shapes
:�
Bgradient_tape/ResNet_entropy_closure/add/add/BroadcastGradientArgsBroadcastGradientArgs;gradient_tape/ResNet_entropy_closure/add/add/Shape:output:0=gradient_tape/ResNet_entropy_closure/add/add/Shape_1:output:0*2
_output_shapes 
:���������:����������
0gradient_tape/ResNet_entropy_closure/add/add/SumSumAddN_4:sum:0Ggradient_tape/ResNet_entropy_closure/add/add/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
4gradient_tape/ResNet_entropy_closure/add/add/ReshapeReshape9gradient_tape/ResNet_entropy_closure/add/add/Sum:output:0;gradient_tape/ResNet_entropy_closure/add/add/Shape:output:0*
T0*(
_output_shapes
:�����������
2gradient_tape/ResNet_entropy_closure/add/add/Sum_1SumAddN_4:sum:0Ggradient_tape/ResNet_entropy_closure/add/add/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
6gradient_tape/ResNet_entropy_closure/add/add/Reshape_1Reshape;gradient_tape/ResNet_entropy_closure/add/add/Sum_1:output:0=gradient_tape/ResNet_entropy_closure/add/add/Shape_1:output:0*
T0*(
_output_shapes
:�����������
Hgradient_tape/ResNet_entropy_closure/block_0_layer_0/BiasAdd/BiasAddGradBiasAddGrad?gradient_tape/ResNet_entropy_closure/add/add/Reshape_1:output:0*
T0*
_output_shapes	
:��
Bgradient_tape/ResNet_entropy_closure/block_0_layer_0/MatMul/MatMulMatMul?gradient_tape/ResNet_entropy_closure/add/add/Reshape_1:output:0DResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������*
transpose_b(�
Dgradient_tape/ResNet_entropy_closure/block_0_layer_0/MatMul/MatMul_1MatMul>ResNet_entropy_closure/tf.math.softplus/Softplus:activations:0?gradient_tape/ResNet_entropy_closure/add/add/Reshape_1:output:0*
T0* 
_output_shapes
:
��*
transpose_a(�
=gradient_tape/ResNet_entropy_closure/tf.math.softplus/SigmoidSigmoid3ResNet_entropy_closure/layer_input/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
9gradient_tape/ResNet_entropy_closure/tf.math.softplus/mulMulLgradient_tape/ResNet_entropy_closure/block_0_layer_0/MatMul/MatMul:product:0Agradient_tape/ResNet_entropy_closure/tf.math.softplus/Sigmoid:y:0*
T0*(
_output_shapes
:�����������
AddN_5AddN=gradient_tape/ResNet_entropy_closure/add/add/Reshape:output:0=gradient_tape/ResNet_entropy_closure/tf.math.softplus/mul:z:0*
N*
T0*(
_output_shapes
:�����������
Dgradient_tape/ResNet_entropy_closure/layer_input/BiasAdd/BiasAddGradBiasAddGradAddN_5:sum:0*
T0*
_output_shapes	
:��
>gradient_tape/ResNet_entropy_closure/layer_input/MatMul/MatMulMatMulAddN_5:sum:0@ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	*
transpose_b(�
@gradient_tape/ResNet_entropy_closure/layer_input/MatMul/MatMul_1MatMul;ResNet_entropy_closure/decorrelation_layer/MatMul:product:0AddN_5:sum:0*
T0*
_output_shapes
:		�*
transpose_a(�
Fgradient_tape/ResNet_entropy_closure/decorrelation_layer/MatMul/MatMulMatMulHgradient_tape/ResNet_entropy_closure/layer_input/MatMul/MatMul:product:0HResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	*
transpose_b(p
?gradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/ShapeShapex*
T0*
_output_shapes
:�
Agradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Shape_1ShapeBResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp:value:0*
T0*
_output_shapes
:�
Ogradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/BroadcastGradientArgsBroadcastGradientArgsHgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Shape:output:0Jgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Shape_1:output:0*2
_output_shapes 
:���������:����������
=gradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/SumSumPgradient_tape/ResNet_entropy_closure/decorrelation_layer/MatMul/MatMul:product:0Tgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
Agradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/ReshapeReshapeFgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Sum:output:0Hgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Shape:output:0*
T0*'
_output_shapes
:���������	�
CastCastJgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Reshape:output:0*

DstT0*

SrcT0*'
_output_shapes
:���������	�
checkedCheckNumericsCast:y:0*
T0*'
_output_shapes
:���������	*d
messageYWinput tensor checking error at alpha = Tensor("Cast:0", shape=(None, 9), dtype=float64)d
checkedandclipped/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped/MinimumMinimumchecked:output:0$checkedandclipped/Minimum/y:output:0*
T0*'
_output_shapes
:���������	\
checkedandclipped/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclippedMaximumcheckedandclipped/Minimum:z:0checkedandclipped/y:output:0*
T0*'
_output_shapes
:���������	d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSliceunknownstrided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:		�*

begin_mask*
end_maskX
Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:X
Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB: T
Tensordot/ShapeShapecheckedandclipped:z:0*
T0*
_output_shapes
:Y
Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2GatherV2Tensordot/Shape:output:0Tensordot/free:output:0 Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2_1GatherV2Tensordot/Shape:output:0Tensordot/axes:output:0"Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:Y
Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: n
Tensordot/ProdProdTensordot/GatherV2:output:0Tensordot/Const:output:0*
T0*
_output_shapes
: [
Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: t
Tensordot/Prod_1ProdTensordot/GatherV2_1:output:0Tensordot/Const_1:output:0*
T0*
_output_shapes
: W
Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concatConcatV2Tensordot/free:output:0Tensordot/axes:output:0Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:y
Tensordot/stackPackTensordot/Prod:output:0Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot/transpose	Transposecheckedandclipped:z:0Tensordot/concat:output:0*
T0*'
_output_shapes
:���������	�
Tensordot/ReshapeReshapeTensordot/transpose:y:0Tensordot/stack:output:0*
T0*0
_output_shapes
:�������������������
Tensordot/MatMulMatMulTensordot/Reshape:output:0strided_slice:output:0*
T0*(
_output_shapes
:����������\
Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�Y
Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concat_1ConcatV2Tensordot/GatherV2:output:0Tensordot/Const_2:output:0 Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
	TensordotReshapeTensordot/MatMul:product:0Tensordot/concat_1:output:0*
T0*(
_output_shapes
:����������Q
ExpExpTensordot:output:0*
T0*(
_output_shapes
:����������Z
Tensordot_1/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_1/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_1/ShapeShapeExp:y:0*
T0*
_output_shapes
:[
Tensordot_1/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2GatherV2Tensordot_1/Shape:output:0Tensordot_1/free:output:0"Tensordot_1/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_1/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2_1GatherV2Tensordot_1/Shape:output:0Tensordot_1/axes:output:0$Tensordot_1/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_1/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_1/ProdProdTensordot_1/GatherV2:output:0Tensordot_1/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_1/Prod_1ProdTensordot_1/GatherV2_1:output:0Tensordot_1/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concatConcatV2Tensordot_1/free:output:0Tensordot_1/axes:output:0 Tensordot_1/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_1/stackPackTensordot_1/Prod:output:0Tensordot_1/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_1/transpose	TransposeExp:y:0Tensordot_1/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_1/ReshapeReshapeTensordot_1/transpose:y:0Tensordot_1/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_1/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_1/transpose_1	Transposetensordot_1_b%Tensordot_1/transpose_1/perm:output:0*
T0*
_output_shapes
:	��
Tensordot_1/MatMulMatMulTensordot_1/Reshape:output:0Tensordot_1/transpose_1:y:0*
T0*'
_output_shapes
:���������]
Tensordot_1/Const_2Const*
_output_shapes
:*
dtype0*
valueB:[
Tensordot_1/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concat_1ConcatV2Tensordot_1/GatherV2:output:0Tensordot_1/Const_2:output:0"Tensordot_1/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_1ReshapeTensordot_1/MatMul:product:0Tensordot_1/concat_1:output:0*
T0*'
_output_shapes
:���������R
LogLogTensordot_1:output:0*
T0*'
_output_shapes
:���������E
NegNegLog:y:0*
T0*'
_output_shapes
:���������M
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :v
concatConcatV2Neg:y:0Cast:y:0concat/axis:output:0*
N*
T0*'
_output_shapes
:���������
�
	checked_1CheckNumericsconcat:output:0^checked*
T0*'
_output_shapes
:���������
*(
messageinput tensor checking errorf
checkedandclipped_1/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped_1/MinimumMinimumchecked_1:output:0&checkedandclipped_1/Minimum/y:output:0*
T0*'
_output_shapes
:���������
^
checkedandclipped_1/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclipped_1Maximumcheckedandclipped_1/Minimum:z:0checkedandclipped_1/y:output:0*
T0*'
_output_shapes
:���������
Z
Tensordot_2/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_2/freeConst*
_output_shapes
:*
dtype0*
valueB: X
Tensordot_2/ShapeShapecheckedandclipped_1:z:0*
T0*
_output_shapes
:[
Tensordot_2/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2GatherV2Tensordot_2/Shape:output:0Tensordot_2/free:output:0"Tensordot_2/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_2/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2_1GatherV2Tensordot_2/Shape:output:0Tensordot_2/axes:output:0$Tensordot_2/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_2/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_2/ProdProdTensordot_2/GatherV2:output:0Tensordot_2/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_2/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_2/Prod_1ProdTensordot_2/GatherV2_1:output:0Tensordot_2/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_2/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concatConcatV2Tensordot_2/free:output:0Tensordot_2/axes:output:0 Tensordot_2/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_2/stackPackTensordot_2/Prod:output:0Tensordot_2/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2/transpose	Transposecheckedandclipped_1:z:0Tensordot_2/concat:output:0*
T0*'
_output_shapes
:���������
�
Tensordot_2/ReshapeReshapeTensordot_2/transpose:y:0Tensordot_2/stack:output:0*
T0*0
_output_shapes
:������������������v
Tensordot_2/MatMulMatMulTensordot_2/Reshape:output:0unknown*
T0*(
_output_shapes
:����������^
Tensordot_2/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�[
Tensordot_2/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concat_1ConcatV2Tensordot_2/GatherV2:output:0Tensordot_2/Const_2:output:0"Tensordot_2/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2ReshapeTensordot_2/MatMul:product:0Tensordot_2/concat_1:output:0*
T0*(
_output_shapes
:����������U
Exp_1ExpTensordot_2:output:0*
T0*(
_output_shapes
:����������W
MulMul	Exp_1:y:0tensordot_1_b*
T0*(
_output_shapes
:����������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSliceunknownstrided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
:	
�*

begin_mask*
end_maskZ
Tensordot_3/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_3/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_3/ShapeShapeMul:z:0*
T0*
_output_shapes
:[
Tensordot_3/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2GatherV2Tensordot_3/Shape:output:0Tensordot_3/free:output:0"Tensordot_3/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_3/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2_1GatherV2Tensordot_3/Shape:output:0Tensordot_3/axes:output:0$Tensordot_3/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_3/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_3/ProdProdTensordot_3/GatherV2:output:0Tensordot_3/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_3/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_3/Prod_1ProdTensordot_3/GatherV2_1:output:0Tensordot_3/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_3/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concatConcatV2Tensordot_3/free:output:0Tensordot_3/axes:output:0 Tensordot_3/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_3/stackPackTensordot_3/Prod:output:0Tensordot_3/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_3/transpose	TransposeMul:z:0Tensordot_3/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_3/ReshapeReshapeTensordot_3/transpose:y:0Tensordot_3/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_3/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_3/transpose_1	Transposestrided_slice_1:output:0%Tensordot_3/transpose_1/perm:output:0*
T0*
_output_shapes
:	�
�
Tensordot_3/MatMulMatMulTensordot_3/Reshape:output:0Tensordot_3/transpose_1:y:0*
T0*'
_output_shapes
:���������
]
Tensordot_3/Const_2Const*
_output_shapes
:*
dtype0*
valueB:
[
Tensordot_3/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concat_1ConcatV2Tensordot_3/GatherV2:output:0Tensordot_3/Const_2:output:0"Tensordot_3/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_3ReshapeTensordot_3/MatMul:product:0Tensordot_3/concat_1:output:0*
T0*'
_output_shapes
:���������
X
Mul_1Mulmul_1_xconcat:output:0*
T0*'
_output_shapes
:���������
_
addAddV2Tensordot_3:output:0	Mul_1:z:0*
T0*'
_output_shapes
:���������
f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSliceadd:z:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������	*

begin_mask*
end_mask�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOpAresnet_entropy_closure_layer_input_matmul_readvariableop_resource*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOpBresnet_entropy_closure_layer_input_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_0_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_0_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_1_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_1_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_2_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_2_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_3_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_3_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_4_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_4_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpEresnet_entropy_closure_block_5_layer_0_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpFresnet_entropy_closure_block_5_layer_0_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpBresnet_entropy_closure_dense_output_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
IdentityIdentity4ResNet_entropy_closure/dense_output/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������

Identity_1IdentityJgradient_tape/ResNet_entropy_closure/mean_shift_layer/sub/Reshape:output:0^NoOp*
T0*'
_output_shapes
:���������	i

Identity_2Identitystrided_slice_2:output:0^NoOp*
T0*'
_output_shapes
:���������	�
NoOpNoOp>^ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp>^ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp=^ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOpA^ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp;^ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp:^ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp:^ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp9^ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp;^ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp^checked
^checked_16^dense_output/kernel/Regularizer/Square/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
2~
=ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_0_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_0_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_1_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_1_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_2_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_2_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_3_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_3_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_4_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_4_layer_0/MatMul/ReadVariableOp2~
=ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp=ResNet_entropy_closure/block_5_layer_0/BiasAdd/ReadVariableOp2|
<ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp<ResNet_entropy_closure/block_5_layer_0/MatMul/ReadVariableOp2�
@ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp@ResNet_entropy_closure/decorrelation_layer/MatMul/ReadVariableOp2x
:ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp:ResNet_entropy_closure/dense_output/BiasAdd/ReadVariableOp2v
9ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp9ResNet_entropy_closure/dense_output/MatMul/ReadVariableOp2v
9ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp9ResNet_entropy_closure/layer_input/BiasAdd/ReadVariableOp2t
8ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp8ResNet_entropy_closure/layer_input/MatMul/ReadVariableOp2x
:ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp:ResNet_entropy_closure/mean_shift_layer/sub/ReadVariableOp2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2
checkedchecked2
	checked_1	checked_12n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp:J F
'
_output_shapes
:���������	

_user_specified_namex:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
�
b__inference___backward_mean_shift_layer_layer_call_and_return_conditional_losses_93976803_93976821
placeholder#
gradients_sub_grad_shape_inputs1
-gradients_sub_grad_shape_1_sub_readvariableop
identity

identity_1^
gradients/grad_ys_0Identityplaceholder*
T0*'
_output_shapes
:���������	g
gradients/sub_grad/ShapeShapegradients_sub_grad_shape_inputs*
T0*
_output_shapes
:w
gradients/sub_grad/Shape_1Shape-gradients_sub_grad_shape_1_sub_readvariableop*
T0*
_output_shapes
:�
(gradients/sub_grad/BroadcastGradientArgsBroadcastGradientArgs!gradients/sub_grad/Shape:output:0#gradients/sub_grad/Shape_1:output:0*2
_output_shapes 
:���������:����������
gradients/sub_grad/SumSumgradients/grad_ys_0:output:0-gradients/sub_grad/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
gradients/sub_grad/ReshapeReshapegradients/sub_grad/Sum:output:0!gradients/sub_grad/Shape:output:0*
T0*'
_output_shapes
:���������	m
gradients/sub_grad/NegNeggradients/grad_ys_0:output:0*
T0*'
_output_shapes
:���������	�
gradients/sub_grad/Sum_1Sumgradients/sub_grad/Neg:y:0-gradients/sub_grad/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
gradients/sub_grad/Reshape_1Reshape!gradients/sub_grad/Sum_1:output:0#gradients/sub_grad/Shape_1:output:0*
T0*
_output_shapes
:	k
IdentityIdentity#gradients/sub_grad/Reshape:output:0*
T0*'
_output_shapes
:���������	b

Identity_1Identity%gradients/sub_grad/Reshape_1:output:0*
T0*
_output_shapes
:	"
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*?
_input_shapes.
,:���������	:���������	:	*g
forward_function_nameNL__forward_mean_shift_layer_layer_call_and_return_conditional_losses_93976820:- )
'
_output_shapes
:���������	:-)
'
_output_shapes
:���������	: 

_output_shapes
:	
�
�
a__inference___backward_block_2_layer_0_layer_call_and_return_conditional_losses_93976641_93976654
placeholder6
2gradients_matmul_grad_matmul_matmul_readvariableop)
%gradients_matmul_grad_matmul_1_inputs
identity

identity_1

identity_2_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������u
"gradients/BiasAdd_grad/BiasAddGradBiasAddGradgradients/grad_ys_0:output:0*
T0*
_output_shapes	
:��
gradients/MatMul_grad/MatMulMatMulgradients/grad_ys_0:output:02gradients_matmul_grad_matmul_matmul_readvariableop*
T0*(
_output_shapes
:����������*
transpose_b(�
gradients/MatMul_grad/MatMul_1MatMul%gradients_matmul_grad_matmul_1_inputsgradients/grad_ys_0:output:0*
T0* 
_output_shapes
:
��*
transpose_a(o
IdentityIdentity&gradients/MatMul_grad/MatMul:product:0*
T0*(
_output_shapes
:����������k

Identity_1Identity(gradients/MatMul_grad/MatMul_1:product:0*
T0* 
_output_shapes
:
��i

Identity_2Identity+gradients/BiasAdd_grad/BiasAddGrad:output:0*
T0*
_output_shapes	
:�"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*G
_input_shapes6
4:����������:
��:����������*f
forward_function_nameMK__forward_block_2_layer_0_layer_call_and_return_conditional_losses_93976653:. *
(
_output_shapes
:����������:&"
 
_output_shapes
:
��:.*
(
_output_shapes
:����������
�
�
__inference_loss_fn_11_93980514N
?block_4_layer_0_bias_regularizer_square_readvariableop_resource:	�
identity��6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOp?block_4_layer_0_bias_regularizer_square_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: f
IdentityIdentity(block_4_layer_0/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 
NoOpNoOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp
�
�
L__forward_mean_shift_layer_layer_call_and_return_conditional_losses_93976820
inputs_0)
sub_readvariableop_resource:	
identity

inputs
sub_readvariableop��sub/ReadVariableOpj
sub/ReadVariableOpReadVariableOpsub_readvariableop_resource*
_output_shapes
:	*
dtype0b
subSubinputs_0sub/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	V
IdentityIdentitysub:z:0^NoOp*
T0*'
_output_shapes
:���������	[
NoOpNoOp^sub/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"0
sub_readvariableopsub/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime*(
_input_shapes
:���������	: *~
backward_function_namedb__inference___backward_mean_shift_layer_layer_call_and_return_conditional_losses_93976803_939768212(
sub/ReadVariableOpsub/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
?__forward_add_layer_call_and_return_conditional_losses_93976731
inputs_0

inputs_1_0
identity

inputs
inputs_1U
addAddV2inputs_0
inputs_1_0*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"
inputsinputs_0"
inputs_1
inputs_1_0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������*q
backward_function_nameWU__inference___backward_add_layer_call_and_return_conditional_losses_93976715_93976732:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
N__inference_mean_shift_layer_layer_call_and_return_conditional_losses_93979964

inputs)
sub_readvariableop_resource:	
identity��sub/ReadVariableOpj
sub/ReadVariableOpReadVariableOpsub_readvariableop_resource*
_output_shapes
:	*
dtype0`
subSubinputssub/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������	V
IdentityIdentitysub:z:0^NoOp*
T0*'
_output_shapes
:���������	[
NoOpNoOp^sub/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*(
_input_shapes
:���������	: 2(
sub/ReadVariableOpsub/ReadVariableOp:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
o
C__inference_add_1_layer_call_and_return_conditional_losses_93980131
inputs_0
inputs_1
identityS
addAddV2inputs_0inputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
�
T
(__inference_add_1_layer_call_fn_93980125
inputs_0
inputs_1
identity�
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_1_layer_call_and_return_conditional_losses_93975379a
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:R N
(
_output_shapes
:����������
"
_user_specified_name
inputs/0:RN
(
_output_shapes
:����������
"
_user_specified_name
inputs/1
�
�
__inference_loss_fn_14_93980547Q
>dense_output_kernel_regularizer_square_readvariableop_resource:	�
identity��5dense_output/kernel/Regularizer/Square/ReadVariableOp�
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOp>dense_output_kernel_regularizer_square_readvariableop_resource*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: e
IdentityIdentity'dense_output/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: ~
NoOpNoOp6^dense_output/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp
�
�
]__inference___backward_layer_input_layer_call_and_return_conditional_losses_93976764_93976777
placeholder6
2gradients_matmul_grad_matmul_matmul_readvariableop)
%gradients_matmul_grad_matmul_1_inputs
identity

identity_1

identity_2_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������u
"gradients/BiasAdd_grad/BiasAddGradBiasAddGradgradients/grad_ys_0:output:0*
T0*
_output_shapes	
:��
gradients/MatMul_grad/MatMulMatMulgradients/grad_ys_0:output:02gradients_matmul_grad_matmul_matmul_readvariableop*
T0*'
_output_shapes
:���������	*
transpose_b(�
gradients/MatMul_grad/MatMul_1MatMul%gradients_matmul_grad_matmul_1_inputsgradients/grad_ys_0:output:0*
T0*
_output_shapes
:		�*
transpose_a(n
IdentityIdentity&gradients/MatMul_grad/MatMul:product:0*
T0*'
_output_shapes
:���������	j

Identity_1Identity(gradients/MatMul_grad/MatMul_1:product:0*
T0*
_output_shapes
:		�i

Identity_2Identity+gradients/BiasAdd_grad/BiasAddGrad:output:0*
T0*
_output_shapes	
:�"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*E
_input_shapes4
2:����������:		�:���������	*b
forward_function_nameIG__forward_layer_input_layer_call_and_return_conditional_losses_93976776:. *
(
_output_shapes
:����������:%!

_output_shapes
:		�:-)
'
_output_shapes
:���������	
��
�
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976014

inputs'
mean_shift_layer_93975865:	.
decorrelation_layer_93975868:		'
layer_input_93975871:		�#
layer_input_93975873:	�,
block_0_layer_0_93975877:
��'
block_0_layer_0_93975879:	�,
block_1_layer_0_93975884:
��'
block_1_layer_0_93975886:	�,
block_2_layer_0_93975891:
��'
block_2_layer_0_93975893:	�,
block_3_layer_0_93975898:
��'
block_3_layer_0_93975900:	�,
block_4_layer_0_93975905:
��'
block_4_layer_0_93975907:	�,
block_5_layer_0_93975912:
��'
block_5_layer_0_93975914:	�(
dense_output_93975918:	�#
dense_output_93975920:
identity��'block_0_layer_0/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_1_layer_0/StatefulPartitionedCall�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_2_layer_0/StatefulPartitionedCall�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_3_layer_0/StatefulPartitionedCall�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_4_layer_0/StatefulPartitionedCall�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_5_layer_0/StatefulPartitionedCall�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�+decorrelation_layer/StatefulPartitionedCall�$dense_output/StatefulPartitionedCall�5dense_output/kernel/Regularizer/Square/ReadVariableOp�#layer_input/StatefulPartitionedCall�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�(mean_shift_layer/StatefulPartitionedCall�
(mean_shift_layer/StatefulPartitionedCallStatefulPartitionedCallinputsmean_shift_layer_93975865*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *W
fRRP
N__inference_mean_shift_layer_layer_call_and_return_conditional_losses_93975264�
+decorrelation_layer/StatefulPartitionedCallStatefulPartitionedCall1mean_shift_layer/StatefulPartitionedCall:output:0decorrelation_layer_93975868*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *Z
fURS
Q__inference_decorrelation_layer_layer_call_and_return_conditional_losses_93975275�
#layer_input/StatefulPartitionedCallStatefulPartitionedCall4decorrelation_layer/StatefulPartitionedCall:output:0layer_input_93975871layer_input_93975873*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *R
fMRK
I__inference_layer_input_layer_call_and_return_conditional_losses_93975301�
tf.math.softplus/SoftplusSoftplus,layer_input/StatefulPartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_0_layer_0/StatefulPartitionedCallStatefulPartitionedCall'tf.math.softplus/Softplus:activations:0block_0_layer_0_93975877block_0_layer_0_93975879*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_0_layer_0_layer_call_and_return_conditional_losses_93975330�
add/PartitionedCallPartitionedCall,layer_input/StatefulPartitionedCall:output:00block_0_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__inference_add_layer_call_and_return_conditional_losses_93975342x
tf.math.softplus_1/SoftplusSoftplusadd/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_1_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_1/Softplus:activations:0block_1_layer_0_93975884block_1_layer_0_93975886*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_1_layer_0_layer_call_and_return_conditional_losses_93975367�
add_1/PartitionedCallPartitionedCalladd/PartitionedCall:output:00block_1_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_1_layer_call_and_return_conditional_losses_93975379z
tf.math.softplus_2/SoftplusSoftplusadd_1/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_2_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_2/Softplus:activations:0block_2_layer_0_93975891block_2_layer_0_93975893*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_2_layer_0_layer_call_and_return_conditional_losses_93975404�
add_2/PartitionedCallPartitionedCalladd_1/PartitionedCall:output:00block_2_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_2_layer_call_and_return_conditional_losses_93975416z
tf.math.softplus_3/SoftplusSoftplusadd_2/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_3_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_3/Softplus:activations:0block_3_layer_0_93975898block_3_layer_0_93975900*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_3_layer_0_layer_call_and_return_conditional_losses_93975441�
add_3/PartitionedCallPartitionedCalladd_2/PartitionedCall:output:00block_3_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_3_layer_call_and_return_conditional_losses_93975453z
tf.math.softplus_4/SoftplusSoftplusadd_3/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_4_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_4/Softplus:activations:0block_4_layer_0_93975905block_4_layer_0_93975907*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_4_layer_0_layer_call_and_return_conditional_losses_93975478�
add_4/PartitionedCallPartitionedCalladd_3/PartitionedCall:output:00block_4_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_4_layer_call_and_return_conditional_losses_93975490z
tf.math.softplus_5/SoftplusSoftplusadd_4/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_5_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_5/Softplus:activations:0block_5_layer_0_93975912block_5_layer_0_93975914*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_5_layer_0_layer_call_and_return_conditional_losses_93975515�
add_5/PartitionedCallPartitionedCalladd_4/PartitionedCall:output:00block_5_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *L
fGRE
C__inference_add_5_layer_call_and_return_conditional_losses_93975527�
$dense_output/StatefulPartitionedCallStatefulPartitionedCalladd_5/PartitionedCall:output:0dense_output_93975918dense_output_93975920*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *S
fNRL
J__inference_dense_output_layer_call_and_return_conditional_losses_93975545�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93975871*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93975873*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93975877* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93975879*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93975884* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93975886*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93975891* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93975893*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93975898* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93975900*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93975905* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93975907*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93975912* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93975914*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_output_93975918*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: |
IdentityIdentity-dense_output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������

NoOpNoOp(^block_0_layer_0/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_1_layer_0/StatefulPartitionedCall7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_2_layer_0/StatefulPartitionedCall7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_3_layer_0/StatefulPartitionedCall7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_4_layer_0/StatefulPartitionedCall7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_5_layer_0/StatefulPartitionedCall7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp,^decorrelation_layer/StatefulPartitionedCall%^dense_output/StatefulPartitionedCall6^dense_output/kernel/Regularizer/Square/ReadVariableOp$^layer_input/StatefulPartitionedCall3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp)^mean_shift_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 2R
'block_0_layer_0/StatefulPartitionedCall'block_0_layer_0/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_1_layer_0/StatefulPartitionedCall'block_1_layer_0/StatefulPartitionedCall2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_2_layer_0/StatefulPartitionedCall'block_2_layer_0/StatefulPartitionedCall2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_3_layer_0/StatefulPartitionedCall'block_3_layer_0/StatefulPartitionedCall2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_4_layer_0/StatefulPartitionedCall'block_4_layer_0/StatefulPartitionedCall2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_5_layer_0/StatefulPartitionedCall'block_5_layer_0/StatefulPartitionedCall2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2Z
+decorrelation_layer/StatefulPartitionedCall+decorrelation_layer/StatefulPartitionedCall2L
$dense_output/StatefulPartitionedCall$dense_output/StatefulPartitionedCall2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2J
#layer_input/StatefulPartitionedCall#layer_input/StatefulPartitionedCall2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp2T
(mean_shift_layer/StatefulPartitionedCall(mean_shift_layer/StatefulPartitionedCall:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
��
�
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93977124
x-
resnet_entropy_closure_93976405:	1
resnet_entropy_closure_93976407:		2
resnet_entropy_closure_93976409:		�.
resnet_entropy_closure_93976411:	�3
resnet_entropy_closure_93976413:
��.
resnet_entropy_closure_93976415:	�3
resnet_entropy_closure_93976417:
��.
resnet_entropy_closure_93976419:	�3
resnet_entropy_closure_93976421:
��.
resnet_entropy_closure_93976423:	�3
resnet_entropy_closure_93976425:
��.
resnet_entropy_closure_93976427:	�3
resnet_entropy_closure_93976429:
��.
resnet_entropy_closure_93976431:	�3
resnet_entropy_closure_93976433:
��.
resnet_entropy_closure_93976435:	�2
resnet_entropy_closure_93976437:	�-
resnet_entropy_closure_93976439:
unknown
tensordot_1_b
mul_1_x
identity

identity_1

identity_2��.ResNet_entropy_closure/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�checked�	checked_1�5dense_output/kernel/Regularizer/Square/ReadVariableOp�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�
.ResNet_entropy_closure/StatefulPartitionedCallStatefulPartitionedCallxresnet_entropy_closure_93976405resnet_entropy_closure_93976407resnet_entropy_closure_93976409resnet_entropy_closure_93976411resnet_entropy_closure_93976413resnet_entropy_closure_93976415resnet_entropy_closure_93976417resnet_entropy_closure_93976419resnet_entropy_closure_93976421resnet_entropy_closure_93976423resnet_entropy_closure_93976425resnet_entropy_closure_93976427resnet_entropy_closure_93976429resnet_entropy_closure_93976431resnet_entropy_closure_93976433resnet_entropy_closure_93976435resnet_entropy_closure_93976437resnet_entropy_closure_93976439*
Tin
2*3
Tout+
)2'*
_collective_manager_ids
 *�
_output_shapes�
�:���������:	�:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:����������:����������:
��:����������:����������:		�:���������	:		:���������	:���������	:	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *[
fVRT
R__forward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976847v
ones_like/ShapeShape7ResNet_entropy_closure/StatefulPartitionedCall:output:0*
T0*
_output_shapes
:T
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?w
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:����������
PartitionedCallPartitionedCallones_like:output:07ResNet_entropy_closure/StatefulPartitionedCall:output:17ResNet_entropy_closure/StatefulPartitionedCall:output:27ResNet_entropy_closure/StatefulPartitionedCall:output:37ResNet_entropy_closure/StatefulPartitionedCall:output:47ResNet_entropy_closure/StatefulPartitionedCall:output:57ResNet_entropy_closure/StatefulPartitionedCall:output:67ResNet_entropy_closure/StatefulPartitionedCall:output:77ResNet_entropy_closure/StatefulPartitionedCall:output:87ResNet_entropy_closure/StatefulPartitionedCall:output:98ResNet_entropy_closure/StatefulPartitionedCall:output:108ResNet_entropy_closure/StatefulPartitionedCall:output:118ResNet_entropy_closure/StatefulPartitionedCall:output:128ResNet_entropy_closure/StatefulPartitionedCall:output:138ResNet_entropy_closure/StatefulPartitionedCall:output:148ResNet_entropy_closure/StatefulPartitionedCall:output:158ResNet_entropy_closure/StatefulPartitionedCall:output:168ResNet_entropy_closure/StatefulPartitionedCall:output:178ResNet_entropy_closure/StatefulPartitionedCall:output:188ResNet_entropy_closure/StatefulPartitionedCall:output:198ResNet_entropy_closure/StatefulPartitionedCall:output:208ResNet_entropy_closure/StatefulPartitionedCall:output:218ResNet_entropy_closure/StatefulPartitionedCall:output:228ResNet_entropy_closure/StatefulPartitionedCall:output:238ResNet_entropy_closure/StatefulPartitionedCall:output:248ResNet_entropy_closure/StatefulPartitionedCall:output:258ResNet_entropy_closure/StatefulPartitionedCall:output:268ResNet_entropy_closure/StatefulPartitionedCall:output:278ResNet_entropy_closure/StatefulPartitionedCall:output:288ResNet_entropy_closure/StatefulPartitionedCall:output:298ResNet_entropy_closure/StatefulPartitionedCall:output:308ResNet_entropy_closure/StatefulPartitionedCall:output:318ResNet_entropy_closure/StatefulPartitionedCall:output:328ResNet_entropy_closure/StatefulPartitionedCall:output:338ResNet_entropy_closure/StatefulPartitionedCall:output:348ResNet_entropy_closure/StatefulPartitionedCall:output:358ResNet_entropy_closure/StatefulPartitionedCall:output:368ResNet_entropy_closure/StatefulPartitionedCall:output:378ResNet_entropy_closure/StatefulPartitionedCall:output:38*2
Tin+
)2'*
Tout
2*
_collective_manager_ids
 *�
_output_shapes�
�:���������	:	:		:		�:�:
��:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�:* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *q
flRj
h__inference___backward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976445_93976848g
CastCastPartitionedCall:output:0*

DstT0*

SrcT0*'
_output_shapes
:���������	�
checkedCheckNumericsCast:y:0*
T0*'
_output_shapes
:���������	*d
messageYWinput tensor checking error at alpha = Tensor("Cast:0", shape=(None, 9), dtype=float64)d
checkedandclipped/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped/MinimumMinimumchecked:output:0$checkedandclipped/Minimum/y:output:0*
T0*'
_output_shapes
:���������	\
checkedandclipped/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclippedMaximumcheckedandclipped/Minimum:z:0checkedandclipped/y:output:0*
T0*'
_output_shapes
:���������	d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSliceunknownstrided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:		�*

begin_mask*
end_maskX
Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:X
Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB: T
Tensordot/ShapeShapecheckedandclipped:z:0*
T0*
_output_shapes
:Y
Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2GatherV2Tensordot/Shape:output:0Tensordot/free:output:0 Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/GatherV2_1GatherV2Tensordot/Shape:output:0Tensordot/axes:output:0"Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:Y
Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: n
Tensordot/ProdProdTensordot/GatherV2:output:0Tensordot/Const:output:0*
T0*
_output_shapes
: [
Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: t
Tensordot/Prod_1ProdTensordot/GatherV2_1:output:0Tensordot/Const_1:output:0*
T0*
_output_shapes
: W
Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concatConcatV2Tensordot/free:output:0Tensordot/axes:output:0Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:y
Tensordot/stackPackTensordot/Prod:output:0Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot/transpose	Transposecheckedandclipped:z:0Tensordot/concat:output:0*
T0*'
_output_shapes
:���������	�
Tensordot/ReshapeReshapeTensordot/transpose:y:0Tensordot/stack:output:0*
T0*0
_output_shapes
:�������������������
Tensordot/MatMulMatMulTensordot/Reshape:output:0strided_slice:output:0*
T0*(
_output_shapes
:����������\
Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�Y
Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot/concat_1ConcatV2Tensordot/GatherV2:output:0Tensordot/Const_2:output:0 Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
	TensordotReshapeTensordot/MatMul:product:0Tensordot/concat_1:output:0*
T0*(
_output_shapes
:����������Q
ExpExpTensordot:output:0*
T0*(
_output_shapes
:����������Z
Tensordot_1/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_1/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_1/ShapeShapeExp:y:0*
T0*
_output_shapes
:[
Tensordot_1/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2GatherV2Tensordot_1/Shape:output:0Tensordot_1/free:output:0"Tensordot_1/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_1/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/GatherV2_1GatherV2Tensordot_1/Shape:output:0Tensordot_1/axes:output:0$Tensordot_1/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_1/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_1/ProdProdTensordot_1/GatherV2:output:0Tensordot_1/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_1/Prod_1ProdTensordot_1/GatherV2_1:output:0Tensordot_1/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concatConcatV2Tensordot_1/free:output:0Tensordot_1/axes:output:0 Tensordot_1/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_1/stackPackTensordot_1/Prod:output:0Tensordot_1/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_1/transpose	TransposeExp:y:0Tensordot_1/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_1/ReshapeReshapeTensordot_1/transpose:y:0Tensordot_1/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_1/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_1/transpose_1	Transposetensordot_1_b%Tensordot_1/transpose_1/perm:output:0*
T0*
_output_shapes
:	��
Tensordot_1/MatMulMatMulTensordot_1/Reshape:output:0Tensordot_1/transpose_1:y:0*
T0*'
_output_shapes
:���������]
Tensordot_1/Const_2Const*
_output_shapes
:*
dtype0*
valueB:[
Tensordot_1/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_1/concat_1ConcatV2Tensordot_1/GatherV2:output:0Tensordot_1/Const_2:output:0"Tensordot_1/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_1ReshapeTensordot_1/MatMul:product:0Tensordot_1/concat_1:output:0*
T0*'
_output_shapes
:���������R
LogLogTensordot_1:output:0*
T0*'
_output_shapes
:���������E
NegNegLog:y:0*
T0*'
_output_shapes
:���������M
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :v
concatConcatV2Neg:y:0Cast:y:0concat/axis:output:0*
N*
T0*'
_output_shapes
:���������
�
	checked_1CheckNumericsconcat:output:0^checked*
T0*'
_output_shapes
:���������
*(
messageinput tensor checking errorf
checkedandclipped_1/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB 2      I@�
checkedandclipped_1/MinimumMinimumchecked_1:output:0&checkedandclipped_1/Minimum/y:output:0*
T0*'
_output_shapes
:���������
^
checkedandclipped_1/yConst*
_output_shapes
: *
dtype0*
valueB 2      I��
checkedandclipped_1Maximumcheckedandclipped_1/Minimum:z:0checkedandclipped_1/y:output:0*
T0*'
_output_shapes
:���������
Z
Tensordot_2/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_2/freeConst*
_output_shapes
:*
dtype0*
valueB: X
Tensordot_2/ShapeShapecheckedandclipped_1:z:0*
T0*
_output_shapes
:[
Tensordot_2/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2GatherV2Tensordot_2/Shape:output:0Tensordot_2/free:output:0"Tensordot_2/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_2/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/GatherV2_1GatherV2Tensordot_2/Shape:output:0Tensordot_2/axes:output:0$Tensordot_2/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_2/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_2/ProdProdTensordot_2/GatherV2:output:0Tensordot_2/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_2/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_2/Prod_1ProdTensordot_2/GatherV2_1:output:0Tensordot_2/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_2/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concatConcatV2Tensordot_2/free:output:0Tensordot_2/axes:output:0 Tensordot_2/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_2/stackPackTensordot_2/Prod:output:0Tensordot_2/Prod_1:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2/transpose	Transposecheckedandclipped_1:z:0Tensordot_2/concat:output:0*
T0*'
_output_shapes
:���������
�
Tensordot_2/ReshapeReshapeTensordot_2/transpose:y:0Tensordot_2/stack:output:0*
T0*0
_output_shapes
:������������������v
Tensordot_2/MatMulMatMulTensordot_2/Reshape:output:0unknown*
T0*(
_output_shapes
:����������^
Tensordot_2/Const_2Const*
_output_shapes
:*
dtype0*
valueB:�[
Tensordot_2/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_2/concat_1ConcatV2Tensordot_2/GatherV2:output:0Tensordot_2/Const_2:output:0"Tensordot_2/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_2ReshapeTensordot_2/MatMul:product:0Tensordot_2/concat_1:output:0*
T0*(
_output_shapes
:����������U
Exp_1ExpTensordot_2:output:0*
T0*(
_output_shapes
:����������W
MulMul	Exp_1:y:0tensordot_1_b*
T0*(
_output_shapes
:����������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSliceunknownstrided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
:	
�*

begin_mask*
end_maskZ
Tensordot_3/axesConst*
_output_shapes
:*
dtype0*
valueB:Z
Tensordot_3/freeConst*
_output_shapes
:*
dtype0*
valueB: H
Tensordot_3/ShapeShapeMul:z:0*
T0*
_output_shapes
:[
Tensordot_3/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2GatherV2Tensordot_3/Shape:output:0Tensordot_3/free:output:0"Tensordot_3/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:]
Tensordot_3/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/GatherV2_1GatherV2Tensordot_3/Shape:output:0Tensordot_3/axes:output:0$Tensordot_3/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:[
Tensordot_3/ConstConst*
_output_shapes
:*
dtype0*
valueB: t
Tensordot_3/ProdProdTensordot_3/GatherV2:output:0Tensordot_3/Const:output:0*
T0*
_output_shapes
: ]
Tensordot_3/Const_1Const*
_output_shapes
:*
dtype0*
valueB: z
Tensordot_3/Prod_1ProdTensordot_3/GatherV2_1:output:0Tensordot_3/Const_1:output:0*
T0*
_output_shapes
: Y
Tensordot_3/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concatConcatV2Tensordot_3/free:output:0Tensordot_3/axes:output:0 Tensordot_3/concat/axis:output:0*
N*
T0*
_output_shapes
:
Tensordot_3/stackPackTensordot_3/Prod:output:0Tensordot_3/Prod_1:output:0*
N*
T0*
_output_shapes
:{
Tensordot_3/transpose	TransposeMul:z:0Tensordot_3/concat:output:0*
T0*(
_output_shapes
:�����������
Tensordot_3/ReshapeReshapeTensordot_3/transpose:y:0Tensordot_3/stack:output:0*
T0*0
_output_shapes
:������������������m
Tensordot_3/transpose_1/permConst*
_output_shapes
:*
dtype0*
valueB"       �
Tensordot_3/transpose_1	Transposestrided_slice_1:output:0%Tensordot_3/transpose_1/perm:output:0*
T0*
_output_shapes
:	�
�
Tensordot_3/MatMulMatMulTensordot_3/Reshape:output:0Tensordot_3/transpose_1:y:0*
T0*'
_output_shapes
:���������
]
Tensordot_3/Const_2Const*
_output_shapes
:*
dtype0*
valueB:
[
Tensordot_3/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : �
Tensordot_3/concat_1ConcatV2Tensordot_3/GatherV2:output:0Tensordot_3/Const_2:output:0"Tensordot_3/concat_1/axis:output:0*
N*
T0*
_output_shapes
:�
Tensordot_3ReshapeTensordot_3/MatMul:product:0Tensordot_3/concat_1:output:0*
T0*'
_output_shapes
:���������
X
Mul_1Mulmul_1_xconcat:output:0*
T0*'
_output_shapes
:���������
_
addAddV2Tensordot_3:output:0	Mul_1:z:0*
T0*'
_output_shapes
:���������
f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSliceadd:z:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������	*

begin_mask*
end_mask�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976409*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976411*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976413* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976415*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976417* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976419*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976421* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976423*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976425* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976427*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976429* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976431*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976433* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976435*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpresnet_entropy_closure_93976437*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
IdentityIdentity7ResNet_entropy_closure/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������i

Identity_1IdentityPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������	i

Identity_2Identitystrided_slice_2:output:0^NoOp*
T0*'
_output_shapes
:���������	�
NoOpNoOp/^ResNet_entropy_closure/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp^checked
^checked_16^dense_output/kernel/Regularizer/Square/ReadVariableOp3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
2`
.ResNet_entropy_closure/StatefulPartitionedCall.ResNet_entropy_closure/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2
checkedchecked2
	checked_1	checked_12n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp:J F
'
_output_shapes
:���������	

_user_specified_namex:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

�
�
&__inference_signature_wrapper_93978582
input_1
unknown:	
	unknown_0:		
	unknown_1:		�
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:
��

unknown_12:	�

unknown_13:
��

unknown_14:	�

unknown_15:	�

unknown_16:

unknown_17

unknown_18

unknown_19
identity

identity_1

identity_2��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19*!
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������	:���������	*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *,
f'R%
#__inference__wrapped_model_93975250o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������	q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������	`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*j
_input_shapesY
W:���������	: : : : : : : : : : : : : : : : : : :	
�:	�:
22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1:%!

_output_shapes
:	
�:%!

_output_shapes
:	�:$ 

_output_shapes

:

��
�
R__forward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93977420

inputs'
mean_shift_layer_93975865:	.
decorrelation_layer_93975868:		'
layer_input_93975871:		�#
layer_input_93975873:	�,
block_0_layer_0_93975877:
��'
block_0_layer_0_93975879:	�,
block_1_layer_0_93975884:
��'
block_1_layer_0_93975886:	�,
block_2_layer_0_93975891:
��'
block_2_layer_0_93975893:	�,
block_3_layer_0_93975898:
��'
block_3_layer_0_93975900:	�,
block_4_layer_0_93975905:
��'
block_4_layer_0_93975907:	�,
block_5_layer_0_93975912:
��'
block_5_layer_0_93975914:	�(
dense_output_93975918:	�#
dense_output_93975920:
identity(
$dense_output_statefulpartitionedcall*
&dense_output_statefulpartitionedcall_0
add_5_partitionedcall
add_5_partitionedcall_0+
'block_5_layer_0_statefulpartitionedcall-
)block_5_layer_0_statefulpartitionedcall_0
add_4_partitionedcall
add_4_partitionedcall_0
add_4_partitionedcall_1+
'block_4_layer_0_statefulpartitionedcall-
)block_4_layer_0_statefulpartitionedcall_0
add_3_partitionedcall
add_3_partitionedcall_0
add_3_partitionedcall_1+
'block_3_layer_0_statefulpartitionedcall-
)block_3_layer_0_statefulpartitionedcall_0
add_2_partitionedcall
add_2_partitionedcall_0
add_2_partitionedcall_1+
'block_2_layer_0_statefulpartitionedcall-
)block_2_layer_0_statefulpartitionedcall_0
add_1_partitionedcall
add_1_partitionedcall_0
add_1_partitionedcall_1+
'block_1_layer_0_statefulpartitionedcall-
)block_1_layer_0_statefulpartitionedcall_0
add_partitionedcall
add_partitionedcall_0
add_partitionedcall_1+
'block_0_layer_0_statefulpartitionedcall-
)block_0_layer_0_statefulpartitionedcall_0'
#layer_input_statefulpartitionedcall)
%layer_input_statefulpartitionedcall_0)
%layer_input_statefulpartitionedcall_1/
+decorrelation_layer_statefulpartitionedcall1
-decorrelation_layer_statefulpartitionedcall_0,
(mean_shift_layer_statefulpartitionedcall.
*mean_shift_layer_statefulpartitionedcall_0��'block_0_layer_0/StatefulPartitionedCall�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_1_layer_0/StatefulPartitionedCall�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_2_layer_0/StatefulPartitionedCall�6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_3_layer_0/StatefulPartitionedCall�6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_4_layer_0/StatefulPartitionedCall�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp�'block_5_layer_0/StatefulPartitionedCall�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp�+decorrelation_layer/StatefulPartitionedCall�$dense_output/StatefulPartitionedCall�5dense_output/kernel/Regularizer/Square/ReadVariableOp�#layer_input/StatefulPartitionedCall�2layer_input/bias/Regularizer/Square/ReadVariableOp�4layer_input/kernel/Regularizer/Square/ReadVariableOp�(mean_shift_layer/StatefulPartitionedCall�
(mean_shift_layer/StatefulPartitionedCallStatefulPartitionedCallinputsmean_shift_layer_93975865*
Tin
2*
Tout
2*
_collective_manager_ids
 *@
_output_shapes.
,:���������	:���������	:	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *U
fPRN
L__forward_mean_shift_layer_layer_call_and_return_conditional_losses_93976820�
+decorrelation_layer/StatefulPartitionedCallStatefulPartitionedCall1mean_shift_layer/StatefulPartitionedCall:output:0decorrelation_layer_93975868*
Tin
2*
Tout
2*
_collective_manager_ids
 *D
_output_shapes2
0:���������	:		:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *X
fSRQ
O__forward_decorrelation_layer_layer_call_and_return_conditional_losses_93976795�
#layer_input/StatefulPartitionedCallStatefulPartitionedCall4decorrelation_layer/StatefulPartitionedCall:output:0layer_input_93975871layer_input_93975873*
Tin
2*
Tout
2*
_collective_manager_ids
 *F
_output_shapes4
2:����������:		�:���������	*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *P
fKRI
G__forward_layer_input_layer_call_and_return_conditional_losses_93976776�
tf.math.softplus/SoftplusSoftplus,layer_input/StatefulPartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_0_layer_0/StatefulPartitionedCallStatefulPartitionedCall'tf.math.softplus/Softplus:activations:0block_0_layer_0_93975877block_0_layer_0_93975879*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_0_layer_0_layer_call_and_return_conditional_losses_93976751�
add/PartitionedCallPartitionedCall,layer_input/StatefulPartitionedCall:output:00block_0_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *H
fCRA
?__forward_add_layer_call_and_return_conditional_losses_93976731x
tf.math.softplus_1/SoftplusSoftplusadd/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_1_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_1/Softplus:activations:0block_1_layer_0_93975884block_1_layer_0_93975886*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_1_layer_0_layer_call_and_return_conditional_losses_93976702�
add_1/PartitionedCallPartitionedCalladd/PartitionedCall:output:00block_1_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_1_layer_call_and_return_conditional_losses_93976682z
tf.math.softplus_2/SoftplusSoftplusadd_1/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_2_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_2/Softplus:activations:0block_2_layer_0_93975891block_2_layer_0_93975893*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_2_layer_0_layer_call_and_return_conditional_losses_93976653�
add_2/PartitionedCallPartitionedCalladd_1/PartitionedCall:output:00block_2_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_2_layer_call_and_return_conditional_losses_93976633z
tf.math.softplus_3/SoftplusSoftplusadd_2/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_3_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_3/Softplus:activations:0block_3_layer_0_93975898block_3_layer_0_93975900*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_3_layer_0_layer_call_and_return_conditional_losses_93976604�
add_3/PartitionedCallPartitionedCalladd_2/PartitionedCall:output:00block_3_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_3_layer_call_and_return_conditional_losses_93976584z
tf.math.softplus_4/SoftplusSoftplusadd_3/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_4_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_4/Softplus:activations:0block_4_layer_0_93975905block_4_layer_0_93975907*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_4_layer_0_layer_call_and_return_conditional_losses_93976555�
add_4/PartitionedCallPartitionedCalladd_3/PartitionedCall:output:00block_4_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_4_layer_call_and_return_conditional_losses_93976535z
tf.math.softplus_5/SoftplusSoftplusadd_4/PartitionedCall:output:0*
T0*(
_output_shapes
:�����������
'block_5_layer_0/StatefulPartitionedCallStatefulPartitionedCall)tf.math.softplus_5/Softplus:activations:0block_5_layer_0_93975912block_5_layer_0_93975914*
Tin
2*
Tout
2*
_collective_manager_ids
 *H
_output_shapes6
4:����������:
��:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *T
fORM
K__forward_block_5_layer_0_layer_call_and_return_conditional_losses_93976506�
add_5/PartitionedCallPartitionedCalladd_4/PartitionedCall:output:00block_5_layer_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *P
_output_shapes>
<:����������:����������:����������* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8� *J
fERC
A__forward_add_5_layer_call_and_return_conditional_losses_93976486�
$dense_output/StatefulPartitionedCallStatefulPartitionedCalladd_5/PartitionedCall:output:0dense_output_93975918dense_output_93975920*
Tin
2*
Tout
2*
_collective_manager_ids
 *F
_output_shapes4
2:���������:	�:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *Q
fLRJ
H__forward_dense_output_layer_call_and_return_conditional_losses_93976461�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93975871*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
2layer_input/bias/Regularizer/Square/ReadVariableOpReadVariableOplayer_input_93975873*
_output_shapes	
:�*
dtype0�
#layer_input/bias/Regularizer/SquareSquare:layer_input/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�l
"layer_input/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
 layer_input/bias/Regularizer/SumSum'layer_input/bias/Regularizer/Square:y:0+layer_input/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: g
"layer_input/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
 layer_input/bias/Regularizer/mulMul+layer_input/bias/Regularizer/mul/x:output:0)layer_input/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93975877* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_0_layer_0_93975879*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93975884* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_1_layer_0_93975886*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93975891* 
_output_shapes
:
��*
dtype0�
)block_2_layer_0/kernel/Regularizer/SquareSquare@block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_2_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_2_layer_0/kernel/Regularizer/SumSum-block_2_layer_0/kernel/Regularizer/Square:y:01block_2_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_2_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_2_layer_0/kernel/Regularizer/mulMul1block_2_layer_0/kernel/Regularizer/mul/x:output:0/block_2_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_2_layer_0_93975893*
_output_shapes	
:�*
dtype0�
'block_2_layer_0/bias/Regularizer/SquareSquare>block_2_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_2_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_2_layer_0/bias/Regularizer/SumSum+block_2_layer_0/bias/Regularizer/Square:y:0/block_2_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_2_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_2_layer_0/bias/Regularizer/mulMul/block_2_layer_0/bias/Regularizer/mul/x:output:0-block_2_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93975898* 
_output_shapes
:
��*
dtype0�
)block_3_layer_0/kernel/Regularizer/SquareSquare@block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_3_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_3_layer_0/kernel/Regularizer/SumSum-block_3_layer_0/kernel/Regularizer/Square:y:01block_3_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_3_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_3_layer_0/kernel/Regularizer/mulMul1block_3_layer_0/kernel/Regularizer/mul/x:output:0/block_3_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_3_layer_0_93975900*
_output_shapes	
:�*
dtype0�
'block_3_layer_0/bias/Regularizer/SquareSquare>block_3_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_3_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_3_layer_0/bias/Regularizer/SumSum+block_3_layer_0/bias/Regularizer/Square:y:0/block_3_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_3_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_3_layer_0/bias/Regularizer/mulMul/block_3_layer_0/bias/Regularizer/mul/x:output:0-block_3_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93975905* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_4_layer_0_93975907*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93975912* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpblock_5_layer_0_93975914*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
5dense_output/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_output_93975918*
_output_shapes
:	�*
dtype0�
&dense_output/kernel/Regularizer/SquareSquare=dense_output/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�v
%dense_output/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
#dense_output/kernel/Regularizer/SumSum*dense_output/kernel/Regularizer/Square:y:0.dense_output/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: j
%dense_output/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
#dense_output/kernel/Regularizer/mulMul.dense_output/kernel/Regularizer/mul/x:output:0,dense_output/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: |
IdentityIdentity-dense_output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������

NoOpNoOp(^block_0_layer_0/StatefulPartitionedCall7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_1_layer_0/StatefulPartitionedCall7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_2_layer_0/StatefulPartitionedCall7^block_2_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_3_layer_0/StatefulPartitionedCall7^block_3_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_4_layer_0/StatefulPartitionedCall7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp(^block_5_layer_0/StatefulPartitionedCall7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp,^decorrelation_layer/StatefulPartitionedCall%^dense_output/StatefulPartitionedCall6^dense_output/kernel/Regularizer/Square/ReadVariableOp$^layer_input/StatefulPartitionedCall3^layer_input/bias/Regularizer/Square/ReadVariableOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp)^mean_shift_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "7
add_1_partitionedcalladd_1/PartitionedCall:output:0"9
add_1_partitionedcall_0add_1/PartitionedCall:output:1"9
add_1_partitionedcall_1add_1/PartitionedCall:output:2"7
add_2_partitionedcalladd_2/PartitionedCall:output:0"9
add_2_partitionedcall_0add_2/PartitionedCall:output:1"9
add_2_partitionedcall_1add_2/PartitionedCall:output:2"7
add_3_partitionedcalladd_3/PartitionedCall:output:0"9
add_3_partitionedcall_0add_3/PartitionedCall:output:1"9
add_3_partitionedcall_1add_3/PartitionedCall:output:2"7
add_4_partitionedcalladd_4/PartitionedCall:output:0"9
add_4_partitionedcall_0add_4/PartitionedCall:output:1"9
add_4_partitionedcall_1add_4/PartitionedCall:output:2"7
add_5_partitionedcalladd_5/PartitionedCall:output:1"9
add_5_partitionedcall_0add_5/PartitionedCall:output:2"3
add_partitionedcalladd/PartitionedCall:output:0"5
add_partitionedcall_0add/PartitionedCall:output:1"5
add_partitionedcall_1add/PartitionedCall:output:2"[
'block_0_layer_0_statefulpartitionedcall0block_0_layer_0/StatefulPartitionedCall:output:1"]
)block_0_layer_0_statefulpartitionedcall_00block_0_layer_0/StatefulPartitionedCall:output:2"[
'block_1_layer_0_statefulpartitionedcall0block_1_layer_0/StatefulPartitionedCall:output:1"]
)block_1_layer_0_statefulpartitionedcall_00block_1_layer_0/StatefulPartitionedCall:output:2"[
'block_2_layer_0_statefulpartitionedcall0block_2_layer_0/StatefulPartitionedCall:output:1"]
)block_2_layer_0_statefulpartitionedcall_00block_2_layer_0/StatefulPartitionedCall:output:2"[
'block_3_layer_0_statefulpartitionedcall0block_3_layer_0/StatefulPartitionedCall:output:1"]
)block_3_layer_0_statefulpartitionedcall_00block_3_layer_0/StatefulPartitionedCall:output:2"[
'block_4_layer_0_statefulpartitionedcall0block_4_layer_0/StatefulPartitionedCall:output:1"]
)block_4_layer_0_statefulpartitionedcall_00block_4_layer_0/StatefulPartitionedCall:output:2"[
'block_5_layer_0_statefulpartitionedcall0block_5_layer_0/StatefulPartitionedCall:output:1"]
)block_5_layer_0_statefulpartitionedcall_00block_5_layer_0/StatefulPartitionedCall:output:2"c
+decorrelation_layer_statefulpartitionedcall4decorrelation_layer/StatefulPartitionedCall:output:1"e
-decorrelation_layer_statefulpartitionedcall_04decorrelation_layer/StatefulPartitionedCall:output:2"U
$dense_output_statefulpartitionedcall-dense_output/StatefulPartitionedCall:output:1"W
&dense_output_statefulpartitionedcall_0-dense_output/StatefulPartitionedCall:output:2"
identityIdentity:output:0"S
#layer_input_statefulpartitionedcall,layer_input/StatefulPartitionedCall:output:0"U
%layer_input_statefulpartitionedcall_0,layer_input/StatefulPartitionedCall:output:1"U
%layer_input_statefulpartitionedcall_1,layer_input/StatefulPartitionedCall:output:2"]
(mean_shift_layer_statefulpartitionedcall1mean_shift_layer/StatefulPartitionedCall:output:1"_
*mean_shift_layer_statefulpartitionedcall_01mean_shift_layer/StatefulPartitionedCall:output:2*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : *�
backward_function_namejh__inference___backward_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93977269_939774212R
'block_0_layer_0/StatefulPartitionedCall'block_0_layer_0/StatefulPartitionedCall2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_1_layer_0/StatefulPartitionedCall'block_1_layer_0/StatefulPartitionedCall2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_2_layer_0/StatefulPartitionedCall'block_2_layer_0/StatefulPartitionedCall2p
6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp6block_2_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_2_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_3_layer_0/StatefulPartitionedCall'block_3_layer_0/StatefulPartitionedCall2p
6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp6block_3_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_3_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_4_layer_0/StatefulPartitionedCall'block_4_layer_0/StatefulPartitionedCall2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp2R
'block_5_layer_0/StatefulPartitionedCall'block_5_layer_0/StatefulPartitionedCall2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp2Z
+decorrelation_layer/StatefulPartitionedCall+decorrelation_layer/StatefulPartitionedCall2L
$dense_output/StatefulPartitionedCall$dense_output/StatefulPartitionedCall2n
5dense_output/kernel/Regularizer/Square/ReadVariableOp5dense_output/kernel/Regularizer/Square/ReadVariableOp2J
#layer_input/StatefulPartitionedCall#layer_input/StatefulPartitionedCall2h
2layer_input/bias/Regularizer/Square/ReadVariableOp2layer_input/bias/Regularizer/Square/ReadVariableOp2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp2T
(mean_shift_layer/StatefulPartitionedCall(mean_shift_layer/StatefulPartitionedCall:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
9__inference_ResNet_entropy_closure_layer_call_fn_93979589

inputs
unknown:	
	unknown_0:		
	unknown_1:		�
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:
��

unknown_12:	�

unknown_13:
��

unknown_14:	�

unknown_15:	�

unknown_16:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *]
fXRV
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93975642o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
K__forward_block_5_layer_0_layer_call_and_return_conditional_losses_93976506
inputs_02
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity
matmul_readvariableop

inputs��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0l
MatMulMatMulinputs_0MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_5_layer_0/kernel/Regularizer/SquareSquare@block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_5_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_5_layer_0/kernel/Regularizer/SumSum-block_5_layer_0/kernel/Regularizer/Square:y:01block_5_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_5_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_5_layer_0/kernel/Regularizer/mulMul1block_5_layer_0/kernel/Regularizer/mul/x:output:0/block_5_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_5_layer_0/bias/Regularizer/SquareSquare>block_5_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_5_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_5_layer_0/bias/Regularizer/SumSum+block_5_layer_0/bias/Regularizer/Square:y:0/block_5_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_5_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_5_layer_0/bias/Regularizer/mulMul/block_5_layer_0/bias/Regularizer/mul/x:output:0-block_5_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_5_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"6
matmul_readvariableopMatMul/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : *}
backward_function_nameca__inference___backward_block_5_layer_0_layer_call_and_return_conditional_losses_93976494_9397650720
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp6block_5_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_5_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
m
C__inference_add_2_layer_call_and_return_conditional_losses_93975416

inputs
inputs_1
identityQ
addAddV2inputsinputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
M__inference_block_4_layer_0_layer_call_and_return_conditional_losses_93975478

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
9__inference_ResNet_entropy_closure_layer_call_fn_93975681
input_1
unknown:	
	unknown_0:		
	unknown_1:		�
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:
��

unknown_12:	�

unknown_13:
��

unknown_14:	�

unknown_15:	�

unknown_16:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8� *]
fXRV
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93975642o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*J
_input_shapes9
7:���������	: : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������	
!
_user_specified_name	input_1
�
�
2__inference_block_1_layer_0_layer_call_fn_93980097

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *V
fQRO
M__inference_block_1_layer_0_layer_call_and_return_conditional_losses_93975367p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
m
C__inference_add_1_layer_call_and_return_conditional_losses_93975379

inputs
inputs_1
identityQ
addAddV2inputsinputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
m
C__inference_add_4_layer_call_and_return_conditional_losses_93975490

inputs
inputs_1
identityQ
addAddV2inputsinputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
3__inference_mean_shift_layer_layer_call_fn_93979957

inputs
unknown:	
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������	*#
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8� *W
fRRP
N__inference_mean_shift_layer_layer_call_and_return_conditional_losses_93975264o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������	`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*(
_input_shapes
:���������	: 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������	
 
_user_specified_nameinputs
�
�
M__inference_block_4_layer_0_layer_call_and_return_conditional_losses_93980284

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_4_layer_0/kernel/Regularizer/SquareSquare@block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_4_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_4_layer_0/kernel/Regularizer/SumSum-block_4_layer_0/kernel/Regularizer/Square:y:01block_4_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_4_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_4_layer_0/kernel/Regularizer/mulMul1block_4_layer_0/kernel/Regularizer/mul/x:output:0/block_4_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_4_layer_0/bias/Regularizer/SquareSquare>block_4_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_4_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_4_layer_0/bias/Regularizer/SumSum+block_4_layer_0/bias/Regularizer/Square:y:0/block_4_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_4_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_4_layer_0/bias/Regularizer/mulMul/block_4_layer_0/bias/Regularizer/mul/x:output:0-block_4_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_4_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp6block_4_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_4_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
K__forward_block_1_layer_0_layer_call_and_return_conditional_losses_93976702
inputs_02
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity
matmul_readvariableop

inputs��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0l
MatMulMatMulinputs_0MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_1_layer_0/bias/Regularizer/SquareSquare>block_1_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_1_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_1_layer_0/bias/Regularizer/SumSum+block_1_layer_0/bias/Regularizer/Square:y:0/block_1_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_1_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_1_layer_0/bias/Regularizer/mulMul/block_1_layer_0/bias/Regularizer/mul/x:output:0-block_1_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_1_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"6
matmul_readvariableopMatMul/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : *}
backward_function_nameca__inference___backward_block_1_layer_0_layer_call_and_return_conditional_losses_93976690_9397670320
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp6block_1_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
K__forward_block_0_layer_0_layer_call_and_return_conditional_losses_93976751
inputs_02
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity
matmul_readvariableop

inputs��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp�8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0l
MatMulMatMulinputs_0MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_0_layer_0/kernel/Regularizer/SquareSquare@block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_0_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_0_layer_0/kernel/Regularizer/SumSum-block_0_layer_0/kernel/Regularizer/Square:y:01block_0_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_0_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_0_layer_0/kernel/Regularizer/mulMul1block_0_layer_0/kernel/Regularizer/mul/x:output:0/block_0_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
'block_0_layer_0/bias/Regularizer/SquareSquare>block_0_layer_0/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�p
&block_0_layer_0/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: �
$block_0_layer_0/bias/Regularizer/SumSum+block_0_layer_0/bias/Regularizer/Square:y:0/block_0_layer_0/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: k
&block_0_layer_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
$block_0_layer_0/bias/Regularizer/mulMul/block_0_layer_0/bias/Regularizer/mul/x:output:0-block_0_layer_0/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: `
IdentityIdentityBiasAdd:output:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp7^block_0_layer_0/bias/Regularizer/Square/ReadVariableOp9^block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"
inputsinputs_0"6
matmul_readvariableopMatMul/ReadVariableOp:value:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : *}
backward_function_nameca__inference___backward_block_0_layer_0_layer_call_and_return_conditional_losses_93976739_9397675220
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2p
6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp6block_0_layer_0/bias/Regularizer/Square/ReadVariableOp2t
8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_0_layer_0/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_0_93980393P
=layer_input_kernel_regularizer_square_readvariableop_resource:		�
identity��4layer_input/kernel/Regularizer/Square/ReadVariableOp�
4layer_input/kernel/Regularizer/Square/ReadVariableOpReadVariableOp=layer_input_kernel_regularizer_square_readvariableop_resource*
_output_shapes
:		�*
dtype0�
%layer_input/kernel/Regularizer/SquareSquare<layer_input/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:		�u
$layer_input/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
"layer_input/kernel/Regularizer/SumSum)layer_input/kernel/Regularizer/Square:y:0-layer_input/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: i
$layer_input/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
"layer_input/kernel/Regularizer/mulMul-layer_input/kernel/Regularizer/mul/x:output:0+layer_input/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: d
IdentityIdentity&layer_input/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: }
NoOpNoOp5^layer_input/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2l
4layer_input/kernel/Regularizer/Square/ReadVariableOp4layer_input/kernel/Regularizer/Square/ReadVariableOp
�
m
C__inference_add_5_layer_call_and_return_conditional_losses_93975527

inputs
inputs_1
identityQ
addAddV2inputsinputs_1*
T0*(
_output_shapes
:����������P
IdentityIdentityadd:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*;
_input_shapes*
(:����������:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs:PL
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_4_93980437U
Ablock_1_layer_0_kernel_regularizer_square_readvariableop_resource:
��
identity��8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp�
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOpReadVariableOpAblock_1_layer_0_kernel_regularizer_square_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
)block_1_layer_0/kernel/Regularizer/SquareSquare@block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��y
(block_1_layer_0/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       �
&block_1_layer_0/kernel/Regularizer/SumSum-block_1_layer_0/kernel/Regularizer/Square:y:01block_1_layer_0/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: m
(block_1_layer_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *��8�
&block_1_layer_0/kernel/Regularizer/mulMul1block_1_layer_0/kernel/Regularizer/mul/x:output:0/block_1_layer_0/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: h
IdentityIdentity*block_1_layer_0/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp9^block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2t
8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp8block_1_layer_0/kernel/Regularizer/Square/ReadVariableOp
�
�
W__inference___backward_add_5_layer_call_and_return_conditional_losses_93976470_93976487
placeholder#
gradients_add_grad_shape_inputs'
#gradients_add_grad_shape_1_inputs_1
identity

identity_1_
gradients/grad_ys_0Identityplaceholder*
T0*(
_output_shapes
:����������g
gradients/add_grad/ShapeShapegradients_add_grad_shape_inputs*
T0*
_output_shapes
:m
gradients/add_grad/Shape_1Shape#gradients_add_grad_shape_1_inputs_1*
T0*
_output_shapes
:�
(gradients/add_grad/BroadcastGradientArgsBroadcastGradientArgs!gradients/add_grad/Shape:output:0#gradients/add_grad/Shape_1:output:0*2
_output_shapes 
:���������:����������
gradients/add_grad/SumSumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r0:0*
T0*
_output_shapes
:�
gradients/add_grad/ReshapeReshapegradients/add_grad/Sum:output:0!gradients/add_grad/Shape:output:0*
T0*(
_output_shapes
:�����������
gradients/add_grad/Sum_1Sumgradients/grad_ys_0:output:0-gradients/add_grad/BroadcastGradientArgs:r1:0*
T0*
_output_shapes
:�
gradients/add_grad/Reshape_1Reshape!gradients/add_grad/Sum_1:output:0#gradients/add_grad/Shape_1:output:0*
T0*(
_output_shapes
:����������l
IdentityIdentity#gradients/add_grad/Reshape:output:0*
T0*(
_output_shapes
:����������p

Identity_1Identity%gradients/add_grad/Reshape_1:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������:����������:����������*\
forward_function_nameCA__forward_add_5_layer_call_and_return_conditional_losses_93976486:. *
(
_output_shapes
:����������:.*
(
_output_shapes
:����������:.*
(
_output_shapes
:����������"�L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
;
input_10
serving_default_input_1:0���������	<
output_10
StatefulPartitionedCall:0���������<
output_20
StatefulPartitionedCall:1���������	<
output_30
StatefulPartitionedCall:2���������	tensorflow/serving/predict:��
�

core_model
	optimizer
loss
	variables
trainable_variables
regularization_losses
	keras_api

signatures
�__call__
+�&call_and_return_all_conditional_losses
�_default_save_signature"
_tf_keras_model
�
	layer-0

layer_with_weights-0

layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer-4
layer_with_weights-3
layer-5
layer-6
layer-7
layer_with_weights-4
layer-8
layer-9
layer-10
layer_with_weights-5
layer-11
layer-12
layer-13
layer_with_weights-6
layer-14
layer-15
layer-16
layer_with_weights-7
layer-17
layer-18
layer-19
layer_with_weights-8
layer-20
layer-21
layer_with_weights-9
layer-22
 	variables
!trainable_variables
"regularization_losses
#	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_network
�
$iter

%beta_1

&beta_2
	'decay
(learning_rate+m�,m�-m�.m�/m�0m�1m�2m�3m�4m�5m�6m�7m�8m�9m�:m�+v�,v�-v�.v�/v�0v�1v�2v�3v�4v�5v�6v�7v�8v�9v�:v�"
	optimizer
 "
trackable_dict_wrapper
�
)0
*1
+2
,3
-4
.5
/6
07
18
29
310
411
512
613
714
815
916
:17"
trackable_list_wrapper
�
+0
,1
-2
.3
/4
05
16
27
38
49
510
611
712
813
914
:15"
trackable_list_wrapper
 "
trackable_list_wrapper
�
;non_trainable_variables

<layers
=metrics
>layer_regularization_losses
?layer_metrics
	variables
trainable_variables
regularization_losses
�__call__
�_default_save_signature
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
-
�serving_default"
signature_map
"
_tf_keras_input_layer
�
)mu
@	variables
Atrainable_variables
Bregularization_losses
C	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�
*
ev_cov_mat
D	variables
Etrainable_variables
Fregularization_losses
G	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�

+kernel
,bias
H	variables
Itrainable_variables
Jregularization_losses
K	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
(
L	keras_api"
_tf_keras_layer
�

-kernel
.bias
M	variables
Ntrainable_variables
Oregularization_losses
P	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�
Q	variables
Rtrainable_variables
Sregularization_losses
T	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
(
U	keras_api"
_tf_keras_layer
�

/kernel
0bias
V	variables
Wtrainable_variables
Xregularization_losses
Y	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�
Z	variables
[trainable_variables
\regularization_losses
]	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
(
^	keras_api"
_tf_keras_layer
�

1kernel
2bias
_	variables
`trainable_variables
aregularization_losses
b	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�
c	variables
dtrainable_variables
eregularization_losses
f	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
(
g	keras_api"
_tf_keras_layer
�

3kernel
4bias
h	variables
itrainable_variables
jregularization_losses
k	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�
l	variables
mtrainable_variables
nregularization_losses
o	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
(
p	keras_api"
_tf_keras_layer
�

5kernel
6bias
q	variables
rtrainable_variables
sregularization_losses
t	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�
u	variables
vtrainable_variables
wregularization_losses
x	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
(
y	keras_api"
_tf_keras_layer
�

7kernel
8bias
z	variables
{trainable_variables
|regularization_losses
}	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�
~	variables
trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�

9kernel
:bias
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�
)0
*1
+2
,3
-4
.5
/6
07
18
29
310
411
512
613
714
815
916
:17"
trackable_list_wrapper
�
+0
,1
-2
.3
/4
05
16
27
38
49
510
611
712
813
914
:15"
trackable_list_wrapper
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13
�14"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
 	variables
!trainable_variables
"regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
:	2Variable
:		2Variable
%:#		�2layer_input/kernel
:�2layer_input/bias
*:(
��2block_0_layer_0/kernel
#:!�2block_0_layer_0/bias
*:(
��2block_1_layer_0/kernel
#:!�2block_1_layer_0/bias
*:(
��2block_2_layer_0/kernel
#:!�2block_2_layer_0/bias
*:(
��2block_3_layer_0/kernel
#:!�2block_3_layer_0/bias
*:(
��2block_4_layer_0/kernel
#:!�2block_4_layer_0/bias
*:(
��2block_5_layer_0/kernel
#:!�2block_5_layer_0/bias
&:$	�2dense_output/kernel
:2dense_output/bias
.
)0
*1"
trackable_list_wrapper
'
0"
trackable_list_wrapper
p
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
'
)0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
@	variables
Atrainable_variables
Bregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
'
*0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
D	variables
Etrainable_variables
Fregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
+0
,1"
trackable_list_wrapper
.
+0
,1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
H	variables
Itrainable_variables
Jregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
"
_generic_user_object
.
-0
.1"
trackable_list_wrapper
.
-0
.1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
M	variables
Ntrainable_variables
Oregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
Q	variables
Rtrainable_variables
Sregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
"
_generic_user_object
.
/0
01"
trackable_list_wrapper
.
/0
01"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
V	variables
Wtrainable_variables
Xregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
Z	variables
[trainable_variables
\regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
"
_generic_user_object
.
10
21"
trackable_list_wrapper
.
10
21"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
_	variables
`trainable_variables
aregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
c	variables
dtrainable_variables
eregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
"
_generic_user_object
.
30
41"
trackable_list_wrapper
.
30
41"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
h	variables
itrainable_variables
jregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
l	variables
mtrainable_variables
nregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
"
_generic_user_object
.
50
61"
trackable_list_wrapper
.
50
61"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
q	variables
rtrainable_variables
sregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
u	variables
vtrainable_variables
wregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
"
_generic_user_object
.
70
81"
trackable_list_wrapper
.
70
81"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
z	variables
{trainable_variables
|regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
~	variables
trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
90
:1"
trackable_list_wrapper
.
90
:1"
trackable_list_wrapper
(
�0"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
)0
*1"
trackable_list_wrapper
�
	0

1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
R

�total

�count
�	variables
�	keras_api"
_tf_keras_metric
R

�total

�count
�	variables
�	keras_api"
_tf_keras_metric
R

�total

�count
�	variables
�	keras_api"
_tf_keras_metric
R

�total

�count
�	variables
�	keras_api"
_tf_keras_metric
c

�total

�count
�
_fn_kwargs
�	variables
�	keras_api"
_tf_keras_metric
c

�total

�count
�
_fn_kwargs
�	variables
�	keras_api"
_tf_keras_metric
c

�total

�count
�
_fn_kwargs
�	variables
�	keras_api"
_tf_keras_metric
c

�total

�count
�
_fn_kwargs
�	variables
�	keras_api"
_tf_keras_metric
c

�total

�count
�
_fn_kwargs
�	variables
�	keras_api"
_tf_keras_metric
c

�total

�count
�
_fn_kwargs
�	variables
�	keras_api"
_tf_keras_metric
'
)0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
'
*0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
(
�0"
trackable_list_wrapper
 "
trackable_dict_wrapper
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
*:(		�2Adam/layer_input/kernel/m
$:"�2Adam/layer_input/bias/m
/:-
��2Adam/block_0_layer_0/kernel/m
(:&�2Adam/block_0_layer_0/bias/m
/:-
��2Adam/block_1_layer_0/kernel/m
(:&�2Adam/block_1_layer_0/bias/m
/:-
��2Adam/block_2_layer_0/kernel/m
(:&�2Adam/block_2_layer_0/bias/m
/:-
��2Adam/block_3_layer_0/kernel/m
(:&�2Adam/block_3_layer_0/bias/m
/:-
��2Adam/block_4_layer_0/kernel/m
(:&�2Adam/block_4_layer_0/bias/m
/:-
��2Adam/block_5_layer_0/kernel/m
(:&�2Adam/block_5_layer_0/bias/m
+:)	�2Adam/dense_output/kernel/m
$:"2Adam/dense_output/bias/m
*:(		�2Adam/layer_input/kernel/v
$:"�2Adam/layer_input/bias/v
/:-
��2Adam/block_0_layer_0/kernel/v
(:&�2Adam/block_0_layer_0/bias/v
/:-
��2Adam/block_1_layer_0/kernel/v
(:&�2Adam/block_1_layer_0/bias/v
/:-
��2Adam/block_2_layer_0/kernel/v
(:&�2Adam/block_2_layer_0/bias/v
/:-
��2Adam/block_3_layer_0/kernel/v
(:&�2Adam/block_3_layer_0/bias/v
/:-
��2Adam/block_4_layer_0/kernel/v
(:&�2Adam/block_4_layer_0/bias/v
/:-
��2Adam/block_5_layer_0/kernel/v
(:&�2Adam/block_5_layer_0/bias/v
+:)	�2Adam/dense_output/kernel/v
$:"2Adam/dense_output/bias/v
�2�
0__inference_sobolev_model_layer_call_fn_93977173
0__inference_sobolev_model_layer_call_fn_93978633
0__inference_sobolev_model_layer_call_fn_93978684
0__inference_sobolev_model_layer_call_fn_93977797�
���
FullArgSpec$
args�
jself
jx

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93979071
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93979458
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93978115
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93978433�
���
FullArgSpec$
args�
jself
jx

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�B�
#__inference__wrapped_model_93975250input_1"�
���
FullArgSpec
args� 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
9__inference_ResNet_entropy_closure_layer_call_fn_93975681
9__inference_ResNet_entropy_closure_layer_call_fn_93979589
9__inference_ResNet_entropy_closure_layer_call_fn_93979630
9__inference_ResNet_entropy_closure_layer_call_fn_93976094�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93979790
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93979950
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976246
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976398�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�B�
&__inference_signature_wrapper_93978582input_1"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
3__inference_mean_shift_layer_layer_call_fn_93979957�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
N__inference_mean_shift_layer_layer_call_and_return_conditional_losses_93979964�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
6__inference_decorrelation_layer_layer_call_fn_93979971�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
Q__inference_decorrelation_layer_layer_call_and_return_conditional_losses_93979978�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
.__inference_layer_input_layer_call_fn_93979999�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
I__inference_layer_input_layer_call_and_return_conditional_losses_93980021�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
2__inference_block_0_layer_0_layer_call_fn_93980042�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
M__inference_block_0_layer_0_layer_call_and_return_conditional_losses_93980064�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
&__inference_add_layer_call_fn_93980070�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
A__inference_add_layer_call_and_return_conditional_losses_93980076�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
2__inference_block_1_layer_0_layer_call_fn_93980097�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
M__inference_block_1_layer_0_layer_call_and_return_conditional_losses_93980119�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
(__inference_add_1_layer_call_fn_93980125�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_add_1_layer_call_and_return_conditional_losses_93980131�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
2__inference_block_2_layer_0_layer_call_fn_93980152�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
M__inference_block_2_layer_0_layer_call_and_return_conditional_losses_93980174�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
(__inference_add_2_layer_call_fn_93980180�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_add_2_layer_call_and_return_conditional_losses_93980186�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
2__inference_block_3_layer_0_layer_call_fn_93980207�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
M__inference_block_3_layer_0_layer_call_and_return_conditional_losses_93980229�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
(__inference_add_3_layer_call_fn_93980235�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_add_3_layer_call_and_return_conditional_losses_93980241�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
2__inference_block_4_layer_0_layer_call_fn_93980262�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
M__inference_block_4_layer_0_layer_call_and_return_conditional_losses_93980284�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
(__inference_add_4_layer_call_fn_93980290�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_add_4_layer_call_and_return_conditional_losses_93980296�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
2__inference_block_5_layer_0_layer_call_fn_93980317�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
M__inference_block_5_layer_0_layer_call_and_return_conditional_losses_93980339�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
(__inference_add_5_layer_call_fn_93980345�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_add_5_layer_call_and_return_conditional_losses_93980351�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
/__inference_dense_output_layer_call_fn_93980366�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
J__inference_dense_output_layer_call_and_return_conditional_losses_93980382�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
__inference_loss_fn_0_93980393�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_1_93980404�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_2_93980415�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_3_93980426�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_4_93980437�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_5_93980448�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_6_93980459�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_7_93980470�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_8_93980481�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_9_93980492�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_10_93980503�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_11_93980514�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_12_93980525�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_13_93980536�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_14_93980547�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
	J
Const
J	
Const_1
J	
Const_2�
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976246u)*+,-./0123456789:8�5
.�+
!�
input_1���������	
p 

 
� "%�"
�
0���������
� �
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93976398u)*+,-./0123456789:8�5
.�+
!�
input_1���������	
p

 
� "%�"
�
0���������
� �
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93979790t)*+,-./0123456789:7�4
-�*
 �
inputs���������	
p 

 
� "%�"
�
0���������
� �
T__inference_ResNet_entropy_closure_layer_call_and_return_conditional_losses_93979950t)*+,-./0123456789:7�4
-�*
 �
inputs���������	
p

 
� "%�"
�
0���������
� �
9__inference_ResNet_entropy_closure_layer_call_fn_93975681h)*+,-./0123456789:8�5
.�+
!�
input_1���������	
p 

 
� "�����������
9__inference_ResNet_entropy_closure_layer_call_fn_93976094h)*+,-./0123456789:8�5
.�+
!�
input_1���������	
p

 
� "�����������
9__inference_ResNet_entropy_closure_layer_call_fn_93979589g)*+,-./0123456789:7�4
-�*
 �
inputs���������	
p 

 
� "�����������
9__inference_ResNet_entropy_closure_layer_call_fn_93979630g)*+,-./0123456789:7�4
-�*
 �
inputs���������	
p

 
� "�����������
#__inference__wrapped_model_93975250�)*+,-./0123456789:���0�-
&�#
!�
input_1���������	
� "���
.
output_1"�
output_1���������
.
output_2"�
output_2���������	
.
output_3"�
output_3���������	�
C__inference_add_1_layer_call_and_return_conditional_losses_93980131�\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "&�#
�
0����������
� �
(__inference_add_1_layer_call_fn_93980125y\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "������������
C__inference_add_2_layer_call_and_return_conditional_losses_93980186�\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "&�#
�
0����������
� �
(__inference_add_2_layer_call_fn_93980180y\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "������������
C__inference_add_3_layer_call_and_return_conditional_losses_93980241�\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "&�#
�
0����������
� �
(__inference_add_3_layer_call_fn_93980235y\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "������������
C__inference_add_4_layer_call_and_return_conditional_losses_93980296�\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "&�#
�
0����������
� �
(__inference_add_4_layer_call_fn_93980290y\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "������������
C__inference_add_5_layer_call_and_return_conditional_losses_93980351�\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "&�#
�
0����������
� �
(__inference_add_5_layer_call_fn_93980345y\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "������������
A__inference_add_layer_call_and_return_conditional_losses_93980076�\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "&�#
�
0����������
� �
&__inference_add_layer_call_fn_93980070y\�Y
R�O
M�J
#� 
inputs/0����������
#� 
inputs/1����������
� "������������
M__inference_block_0_layer_0_layer_call_and_return_conditional_losses_93980064^-.0�-
&�#
!�
inputs����������
� "&�#
�
0����������
� �
2__inference_block_0_layer_0_layer_call_fn_93980042Q-.0�-
&�#
!�
inputs����������
� "������������
M__inference_block_1_layer_0_layer_call_and_return_conditional_losses_93980119^/00�-
&�#
!�
inputs����������
� "&�#
�
0����������
� �
2__inference_block_1_layer_0_layer_call_fn_93980097Q/00�-
&�#
!�
inputs����������
� "������������
M__inference_block_2_layer_0_layer_call_and_return_conditional_losses_93980174^120�-
&�#
!�
inputs����������
� "&�#
�
0����������
� �
2__inference_block_2_layer_0_layer_call_fn_93980152Q120�-
&�#
!�
inputs����������
� "������������
M__inference_block_3_layer_0_layer_call_and_return_conditional_losses_93980229^340�-
&�#
!�
inputs����������
� "&�#
�
0����������
� �
2__inference_block_3_layer_0_layer_call_fn_93980207Q340�-
&�#
!�
inputs����������
� "������������
M__inference_block_4_layer_0_layer_call_and_return_conditional_losses_93980284^560�-
&�#
!�
inputs����������
� "&�#
�
0����������
� �
2__inference_block_4_layer_0_layer_call_fn_93980262Q560�-
&�#
!�
inputs����������
� "������������
M__inference_block_5_layer_0_layer_call_and_return_conditional_losses_93980339^780�-
&�#
!�
inputs����������
� "&�#
�
0����������
� �
2__inference_block_5_layer_0_layer_call_fn_93980317Q780�-
&�#
!�
inputs����������
� "������������
Q__inference_decorrelation_layer_layer_call_and_return_conditional_losses_93979978[*/�,
%�"
 �
inputs���������	
� "%�"
�
0���������	
� �
6__inference_decorrelation_layer_layer_call_fn_93979971N*/�,
%�"
 �
inputs���������	
� "����������	�
J__inference_dense_output_layer_call_and_return_conditional_losses_93980382]9:0�-
&�#
!�
inputs����������
� "%�"
�
0���������
� �
/__inference_dense_output_layer_call_fn_93980366P9:0�-
&�#
!�
inputs����������
� "�����������
I__inference_layer_input_layer_call_and_return_conditional_losses_93980021]+,/�,
%�"
 �
inputs���������	
� "&�#
�
0����������
� �
.__inference_layer_input_layer_call_fn_93979999P+,/�,
%�"
 �
inputs���������	
� "�����������=
__inference_loss_fn_0_93980393+�

� 
� "� >
__inference_loss_fn_10_939805035�

� 
� "� >
__inference_loss_fn_11_939805146�

� 
� "� >
__inference_loss_fn_12_939805257�

� 
� "� >
__inference_loss_fn_13_939805368�

� 
� "� >
__inference_loss_fn_14_939805479�

� 
� "� =
__inference_loss_fn_1_93980404,�

� 
� "� =
__inference_loss_fn_2_93980415-�

� 
� "� =
__inference_loss_fn_3_93980426.�

� 
� "� =
__inference_loss_fn_4_93980437/�

� 
� "� =
__inference_loss_fn_5_939804480�

� 
� "� =
__inference_loss_fn_6_939804591�

� 
� "� =
__inference_loss_fn_7_939804702�

� 
� "� =
__inference_loss_fn_8_939804813�

� 
� "� =
__inference_loss_fn_9_939804924�

� 
� "� �
N__inference_mean_shift_layer_layer_call_and_return_conditional_losses_93979964[)/�,
%�"
 �
inputs���������	
� "%�"
�
0���������	
� �
3__inference_mean_shift_layer_layer_call_fn_93979957N)/�,
%�"
 �
inputs���������	
� "����������	�
&__inference_signature_wrapper_93978582�)*+,-./0123456789:���;�8
� 
1�.
,
input_1!�
input_1���������	"���
.
output_1"�
output_1���������
.
output_2"�
output_2���������	
.
output_3"�
output_3���������	�
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93978115�)*+,-./0123456789:���4�1
*�'
!�
input_1���������	
p 
� "j�g
`�]
�
0/0���������
�
0/1���������	
�
0/2���������	
� �
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93978433�)*+,-./0123456789:���4�1
*�'
!�
input_1���������	
p
� "j�g
`�]
�
0/0���������
�
0/1���������	
�
0/2���������	
� �
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93979071�)*+,-./0123456789:���.�+
$�!
�
x���������	
p 
� "j�g
`�]
�
0/0���������
�
0/1���������	
�
0/2���������	
� �
K__inference_sobolev_model_layer_call_and_return_conditional_losses_93979458�)*+,-./0123456789:���.�+
$�!
�
x���������	
p
� "j�g
`�]
�
0/0���������
�
0/1���������	
�
0/2���������	
� �
0__inference_sobolev_model_layer_call_fn_93977173�)*+,-./0123456789:���4�1
*�'
!�
input_1���������	
p 
� "Z�W
�
0���������
�
1���������	
�
2���������	�
0__inference_sobolev_model_layer_call_fn_93977797�)*+,-./0123456789:���4�1
*�'
!�
input_1���������	
p
� "Z�W
�
0���������
�
1���������	
�
2���������	�
0__inference_sobolev_model_layer_call_fn_93978633�)*+,-./0123456789:���.�+
$�!
�
x���������	
p 
� "Z�W
�
0���������
�
1���������	
�
2���������	�
0__inference_sobolev_model_layer_call_fn_93978684�)*+,-./0123456789:���.�+
$�!
�
x���������	
p
� "Z�W
�
0���������
�
1���������	
�
2���������	