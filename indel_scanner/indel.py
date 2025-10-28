from dataclasses import dataclass
from enum import StrEnum
from abc import ABC, abstractmethod
from typing import List, Optional

class INDEL_TYPE(StrEnum):
	INSERTION = 'ins'
	DELETION = 'del'

@dataclass
class Indel(ABC):
	contig:str
	ref_position:int
	length:int
	prefix_context:str
	suffix_context:str
	read_name:str
	type:INDEL_TYPE
	

	@property
	@abstractmethod    
	def tsv_sequence_field(self) -> str:
		pass

	def to_tsv_row(self) -> str:
		"""Defines the common formatting logic for the entire TSV row."""
		return '\t'.join([
			self.contig,
			str(self.ref_position),
			self.type.value,
			str(self.length),
			self.prefix_context,
			self.tsv_sequence_field,#  abstract property 
			self.suffix_context,
			self.read_name
		])

@dataclass
class Insertion(Indel):
	inserted_seq:str

	@property
	def tsv_sequence_field(self) -> str:	
		return self.inserted_seq

	insertion_quality:Optional[List[int]] = None
	prefix_quality:Optional[List[int]] = None
	suffix_quality:Optional[List[int]] = None

	

@dataclass
class Deletion(Indel):
	reference_seq:str
	
	@property
	def tsv_sequence_field(self) -> str:	
		return self.reference_seq
	
	prefix_quality:Optional[List[int]] = None
	suffix_quality:Optional[List[int]] = None
	